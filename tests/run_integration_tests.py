#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Olga Goulko and David Llamas

import contextlib
import math
import pathlib
import subprocess
import sys
import tempfile


STANDARD_OUTPUTS = {
    "result.dat", "energies.dat", "temperatures.dat", "config.dat",
    "initconf.dat", "finconf.dat", "bestconf.dat",
}


class IntegrationTestFailure(RuntimeError):
    pass


class IntegrationTestSuite:
    def __init__(self):
        self.test_count = 0
        self.failures = []

    @contextlib.contextmanager
    def case(self, name):
        self.test_count += 1
        try:
            yield
        except Exception as error:
            self.failures.append((name, error))

    def finish(self):
        if self.failures:
            print(f"integration tests: FAIL "
                  f"({len(self.failures)} of {self.test_count} tests failed)",
                  file=sys.stderr)
            for name, error in self.failures:
                print(f"- {name}: {error}", file=sys.stderr)
            return 1

        print(f"integration tests: PASS ({self.test_count} tests)")
        return 0


def require(condition, message):
    if not condition:
        raise IntegrationTestFailure(message)


def read_result_values(path):
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        label, separator, value = line.partition(": ")
        if separator:
            values[label] = value
    return values


# Independent reference implementation of the delta-function discretizations.
def reference_contribution_energy(normalized_offset, delta_option):
    if normalized_offset >= 2.0:
        return 0.0

    if delta_option == 0:
        return (1.0 + math.cos(math.pi * normalized_offset / 2.0)) / 4.0

    if normalized_offset < 1.0:
        return (17.0 / 48.0 + math.sqrt(3.0) * math.pi / 108.0
                + normalized_offset / 4.0 - normalized_offset ** 2 / 4.0
                + (1.0 - 2.0 * normalized_offset)
                * math.sqrt(1.0 + 12.0 * normalized_offset
                            * (1.0 - normalized_offset)) / 16.0
                - math.sqrt(3.0)
                * math.asin(math.sqrt(3.0)
                            * (2.0 * normalized_offset - 1.0) / 2.0)
                / 12.0)

    return (55.0 / 48.0 - math.sqrt(3.0) * math.pi / 108.0
            - 13.0 * normalized_offset / 12.0 + normalized_offset ** 2 / 4.0
            + (2.0 * normalized_offset - 3.0)
            * math.sqrt(36.0 * normalized_offset - 23.0
                        - 12.0 * normalized_offset ** 2) / 48.0
            + math.sqrt(3.0)
            * math.asin(math.sqrt(3.0) * (2.0 * normalized_offset - 3.0) / 2.0)
            / 36.0)


def direct_pairwise_energy(coordinates, grid_size, hopping_distance, delta_option):
    cell_size = len(coordinates) ** (-1.0 / 3.0)
    energy = 0.0

    for first_index, first_coordinate in enumerate(coordinates):
        first_x = first_coordinate % grid_size
        first_y = (first_coordinate // grid_size) % grid_size
        first_z = first_coordinate // grid_size // grid_size

        for second_coordinate in coordinates[first_index + 1:]:
            second_x = second_coordinate % grid_size
            second_y = (second_coordinate // grid_size) % grid_size
            second_z = second_coordinate // grid_size // grid_size
            distance = cell_size * math.sqrt(
                (first_x - second_x) ** 2
                + (first_y - second_y) ** 2
                + (first_z - second_z) ** 2)
            normalized_offset = abs(hopping_distance - distance) / cell_size
            energy += reference_contribution_energy(normalized_offset, delta_option)

    return energy


def require_full_precision(actual, expected, description):
    require(math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-12),
            f"{description}: {actual:.17g} != {expected:.17g}")


def normalized_probability(energy, hopping_distance, total_spins):
    return energy / (
        2.0 * math.pi * hopping_distance ** 2 * total_spins ** (5.0 / 3.0))


def centered_block_coordinates(grid_size):
    x_start = (grid_size - 4) // 2
    y_start = (grid_size - 5) // 2
    z_start = (grid_size - 5) // 2
    return [
        z * grid_size * grid_size + y * grid_size + x
        for z in range(z_start, z_start + 5)
        for y in range(y_start, y_start + 5)
        for x in range(x_start, x_start + 4)
    ]


def simulation_command(executable, total_spins, grid_size, hopping_distance,
                       initialization, delta_option=0):
    return [
        str(executable),
        "-N", str(total_spins),
        "-gridsize", str(grid_size),
        "-d", str(hopping_distance),
        "-hours", "0",
        "-steps", "1",
        "-tempsteps", "10",
        "-inittemp", "1",
        "-fintemp", "0.1",
        "-annealsteps", "100",
        "-initconf", initialization,
        "-delta", str(delta_option),
        "-randomseed", "12345",
    ]


def run_probability_case(executable, grid_size, hopping_distance, delta_option):
    total_spins = 100
    cell_size = total_spins ** (-1.0 / 3.0)
    template_reach = (grid_size - 1) // 2
    required_reach = math.ceil(hopping_distance / cell_size) + 1
    require(template_reach == required_reach,
            "probability fixture is not at minimum valid template reach")

    coordinates = centered_block_coordinates(grid_size)
    require(len(coordinates) == total_spins,
            "integration fixture has the wrong number of coordinates")

    with tempfile.TemporaryDirectory(
            prefix=f"grasshopper3d-delta{delta_option}-") as directory:
        working_directory = pathlib.Path(directory)
        (working_directory / "initconf.dat").write_text(
            "".join(f"{coordinate}\n" for coordinate in coordinates),
            encoding="utf-8")

        command = simulation_command(
            executable, total_spins, grid_size, hopping_distance,
            "load", delta_option)
        completed = subprocess.run(
            command,
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode == 0,
                f"grasshopper exited with status {completed.returncode}: "
                f"{completed.stderr}")

        result_path = working_directory / "result.dat"
        result_text = result_path.read_text(encoding="utf-8")
        result_values = read_result_values(result_path)
        require("Finished after 0 steps" in result_text,
                "probability test did not remain at the loaded initial state")

        reference_energy = direct_pairwise_energy(
            coordinates, grid_size, hopping_distance, delta_option)
        # Each pair occurs once in reference_energy. The independently derived
        # three-dimensional normalization is P = E/(2*pi*d^2*N^(5/3)).
        reference_probability = normalized_probability(
            reference_energy, hopping_distance, total_spins)

        initial_energy = float(
            (working_directory / "energies.dat")
            .read_text(encoding="utf-8").splitlines()[0])
        configuration_values = (
            (working_directory / "config.dat")
            .read_text(encoding="utf-8").splitlines()[0].split())
        require(len(configuration_values) == total_spins + 1,
                "initial config.dat row has the wrong number of fields")
        configuration_energy = float(configuration_values[-1])
        reported_energy = float(result_values["final energy"])
        reported_probability = float(result_values["final probability"])
        reported_acceptance = float(result_values["Average acceptance ratio"])

        require_full_precision(
            initial_energy, reference_energy, "initial reported energy")
        require_full_precision(
            configuration_energy, reference_energy,
            "raw energy in initial config.dat row")
        require_full_precision(
            reported_energy, reference_energy, "zero-step result energy")
        require_full_precision(
            reported_probability, reference_probability,
            "zero-step result probability")
        require(math.isfinite(reported_acceptance) and reported_acceptance == 0.0,
                "zero-step acceptance ratio was not reported as zero")


def run_rejection_case(executable, total_spins, grid_size, hopping_distance,
                       expected_error, case_prefix):
    with tempfile.TemporaryDirectory(prefix=case_prefix) as directory:
        completed = subprocess.run(
            simulation_command(
                executable, total_spins, grid_size, hopping_distance, "random"),
            cwd=directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode != 0, "invalid geometry was not rejected")
        require(expected_error in completed.stderr,
                f"expected error {expected_error!r}; stderr was {completed.stderr!r}")


def run_malformed_cli_preserves_result(executable):
    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-malformed-cli-") as directory:
        working_directory = pathlib.Path(directory)
        result_path = working_directory / "result.dat"
        original_contents = "existing research result\n"
        result_path.write_text(original_contents, encoding="utf-8")

        completed = subprocess.run(
            [str(executable), "-d", "not-a-number"],
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode != 0, "malformed CLI input was accepted")
        require(result_path.read_text(encoding="utf-8") == original_contents,
                "malformed CLI input changed existing result.dat")


def small_annealing_command(executable, initialization, temperature_round_steps,
                            annealing_steps):
    return [
        str(executable),
        "-N", "2",
        "-gridsize", "5",
        "-d", "0.1",
        "-hours", "1",
        "-steps", "1",
        "-tempsteps", str(temperature_round_steps),
        "-inittemp", "1",
        "-fintemp", "0.1",
        "-annealsteps", str(annealing_steps),
        "-initconf", initialization,
        "-delta", "0",
        "-randomseed", "12345",
    ]


def run_short_annealing_schedule(executable):
    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-short-annealing-") as directory:
        completed = subprocess.run(
            small_annealing_command(executable, "random", 1, 1),
            cwd=directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode == 0,
                f"short annealing schedule failed: {completed.stderr}")
        result_text = (pathlib.Path(directory) / "result.dat").read_text(
            encoding="utf-8")
        require("Finished after 1 steps" in result_text,
                "short annealing schedule did not complete its proposal")


def run_partial_round_acceptance_case(executable):
    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-partial-acceptance-") as directory:
        working_directory = pathlib.Path(directory)
        initial_coordinates = [0, 1]
        (working_directory / "initconf.dat").write_text(
            "".join(f"{coordinate}\n" for coordinate in initial_coordinates),
            encoding="utf-8")

        completed = subprocess.run(
            small_annealing_command(executable, "load", 10, 100),
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode == 0,
                f"partial-round run failed: {completed.stderr}")

        result_path = working_directory / "result.dat"
        result_values = read_result_values(result_path)
        require("Finished after 1 steps" in result_path.read_text(encoding="utf-8"),
                "partial-round run did not stop after one proposal")
        final_coordinates = [
            int(value) for value in
            (working_directory / "finconf.dat")
            .read_text(encoding="utf-8").split()
        ]
        accepted_moves = int(set(final_coordinates) != set(initial_coordinates))
        reported_acceptance = float(result_values["Average acceptance ratio"])
        require(reported_acceptance == float(accepted_moves),
                "partial-round acceptance ratio disagrees with configuration change")


def output_policy_command(executable, initialization="random", overwrite=None):
    command = [
        str(executable),
        "-N", "2",
        "-gridsize", "5",
        "-d", "0.1",
        "-hours", "0",
        "-steps", "1",
        "-tempsteps", "10",
        "-inittemp", "1",
        "-fintemp", "0.1",
        "-annealsteps", "100",
        "-initconf", initialization,
        "-delta", "0",
        "-randomseed", "12345",
    ]
    if overwrite is not None:
        command.extend(["-overwrite", str(overwrite)])
    return command


def run_default_output_rejection(executable):
    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-no-overwrite-") as directory:
        working_directory = pathlib.Path(directory)
        stale_config = working_directory / "config.dat"
        sentinel = "configuration from an earlier run\n"
        stale_config.write_text(sentinel, encoding="utf-8")

        completed = subprocess.run(
            output_policy_command(executable),
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode != 0,
                "default overwrite policy accepted an existing output")
        require("Output artifact already exists: config.dat" in completed.stderr,
                "default overwrite rejection did not identify config.dat")
        require(stale_config.read_text(encoding="utf-8") == sentinel,
                "default overwrite rejection changed config.dat")
        require(not (working_directory / "result.dat").exists(),
                "default overwrite rejection created result.dat")


def run_explicit_overwrite(executable):
    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-overwrite-") as directory:
        working_directory = pathlib.Path(directory)
        stale_result = working_directory / "result.dat"
        stale_config = working_directory / "config.dat"
        stale_result.write_text("result from an earlier run\n", encoding="utf-8")
        stale_config.write_text("stale configuration output\n", encoding="utf-8")

        completed = subprocess.run(
            output_policy_command(executable, overwrite=1),
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode == 0,
                f"overwrite run failed: {completed.stderr}")
        require("3D Grasshopper with Simulated Annealing" in
                stale_result.read_text(encoding="utf-8"),
                "overwrite run did not replace result.dat")
        require("stale configuration output" not in
                stale_config.read_text(encoding="utf-8"),
                "overwrite run did not replace config.dat")

        require(all((working_directory / filename).is_file()
                    for filename in STANDARD_OUTPUTS),
                "overwrite run did not create the complete standard output set")


def run_load_overwrite_preserves_input(executable):
    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-load-overwrite-") as directory:
        working_directory = pathlib.Path(directory)
        initial_configuration = working_directory / "initconf.dat"
        initial_contents = "0\n1\n"
        initial_configuration.write_text(initial_contents, encoding="utf-8")
        (working_directory / "config.dat").write_text(
            "stale configuration output\n", encoding="utf-8")

        completed = subprocess.run(
            output_policy_command(executable, initialization="load", overwrite=1),
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode == 0,
                f"load overwrite run failed: {completed.stderr}")
        require(initial_configuration.read_text(encoding="utf-8") == initial_contents,
                "load overwrite run changed initconf.dat")
        require("stale configuration output" not in
                (working_directory / "config.dat").read_text(encoding="utf-8"),
                "load overwrite run did not replace config.dat")


def run_output_cleanup_preflight_failure(executable):
    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-cleanup-failure-") as directory:
        working_directory = pathlib.Path(directory)
        stale_result = working_directory / "result.dat"
        sentinel = "result from an earlier run\n"
        stale_result.write_text(sentinel, encoding="utf-8")
        blocked_artifact = working_directory / "config.dat"
        blocked_artifact.mkdir()

        completed = subprocess.run(
            output_policy_command(executable, overwrite=1),
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode != 0,
                "output cleanup preflight failure was ignored")
        require("Output artifact is a directory: config.dat" in completed.stderr,
                "cleanup preflight did not identify directory artifact")
        require(stale_result.read_text(encoding="utf-8") == sentinel,
                "cleanup preflight changed an earlier output")
        require(blocked_artifact.is_dir(),
                "cleanup preflight removed the directory artifact")


def run_invalid_load_preserves_existing_output(executable):
    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-invalid-load-sentinel-") as directory:
        working_directory = pathlib.Path(directory)
        (working_directory / "initconf.dat").write_text("0\n0\n", encoding="utf-8")
        result_path = working_directory / "result.dat"
        sentinel = "result from an earlier run\n"
        result_path.write_text(sentinel, encoding="utf-8")

        completed = subprocess.run(
            output_policy_command(executable, initialization="load", overwrite=1),
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode != 0,
                "duplicate loaded coordinate was accepted")
        require("Duplicate coordinate" in completed.stderr,
                "invalid loaded configuration did not report its duplicate")
        require(result_path.read_text(encoding="utf-8") == sentinel,
                "invalid loaded configuration changed existing result.dat")


def run_invalid_load_creates_no_outputs(executable):
    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-invalid-load-clean-") as directory:
        working_directory = pathlib.Path(directory)
        (working_directory / "initconf.dat").write_text("0\n", encoding="utf-8")

        completed = subprocess.run(
            output_policy_command(executable, initialization="load", overwrite=1),
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode != 0,
                "insufficient loaded configuration was accepted")
        require("Invalid or insufficient data" in completed.stderr,
                "invalid loaded configuration did not report insufficient data")
        created_outputs = [
            filename for filename in STANDARD_OUTPUTS
            if filename != "initconf.dat" and (working_directory / filename).exists()
        ]
        require(not created_outputs,
                f"invalid loaded configuration created outputs: {created_outputs}")


def read_configuration(path, expected_count, grid_volume):
    coordinates = [int(value) for value in path.read_text(encoding="utf-8").split()]
    require(len(coordinates) == expected_count,
            f"{path.name} contains {len(coordinates)} coordinates, expected {expected_count}")
    require(len(set(coordinates)) == expected_count,
            f"{path.name} contains duplicate coordinates")
    require(all(0 <= coordinate < grid_volume for coordinate in coordinates),
            f"{path.name} contains an out-of-range coordinate")
    return coordinates


def run_final_state_energy_audit(executable):
    total_spins = 6
    grid_size = 9
    grid_volume = grid_size ** 3
    hopping_distance = 1.2
    delta_option = 0
    initial_coordinates = [355, 363, 364, 365, 373, 445]

    with tempfile.TemporaryDirectory(
            prefix="grasshopper3d-final-energy-audit-") as directory:
        working_directory = pathlib.Path(directory)
        (working_directory / "initconf.dat").write_text(
            "".join(f"{coordinate}\n" for coordinate in initial_coordinates),
            encoding="utf-8")
        command = [
            str(executable),
            "-N", str(total_spins),
            "-gridsize", str(grid_size),
            "-d", str(hopping_distance),
            "-hours", "1",
            "-steps", "40",
            "-tempsteps", "100",
            "-inittemp", "10",
            "-fintemp", "0.1",
            "-annealsteps", "100",
            "-initconf", "load",
            "-delta", str(delta_option),
            "-randomseed", "24680",
        ]
        completed = subprocess.run(
            command,
            cwd=working_directory,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        require(completed.returncode == 0,
                f"final-state energy audit run failed: {completed.stderr}")

        result_values = read_result_values(working_directory / "result.dat")
        final_coordinates = read_configuration(
            working_directory / "finconf.dat", total_spins, grid_volume)
        best_coordinates = read_configuration(
            working_directory / "bestconf.dat", total_spins, grid_volume)
        require(set(final_coordinates) != set(initial_coordinates),
                "energy audit fixture did not produce an accepted move")

        final_energy = direct_pairwise_energy(
            final_coordinates, grid_size, hopping_distance, delta_option)
        best_energy = direct_pairwise_energy(
            best_coordinates, grid_size, hopping_distance, delta_option)
        final_probability = normalized_probability(
            final_energy, hopping_distance, total_spins)
        best_probability = normalized_probability(
            best_energy, hopping_distance, total_spins)

        require_full_precision(
            float(result_values["final energy"]), final_energy,
            "audited final energy")
        require_full_precision(
            float(result_values["best energy"]), best_energy,
            "audited best energy")
        require_full_precision(
            float(result_values["final probability"]), final_probability,
            "audited final probability")
        require_full_precision(
            float(result_values["best probability"]), best_probability,
            "audited best probability")


def main():
    if len(sys.argv) != 2:
        print(f"usage: {pathlib.Path(sys.argv[0]).name} GRASSHOPPER_EXECUTABLE",
              file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    suite = IntegrationTestSuite()

    try:
        require(executable.is_file(), f"executable not found: {executable}")
        reach_cases = (
            ("minimum even template reach", 10, 0.5),
            ("minimum odd template reach", 11, 0.75),
        )
        for reach_name, grid_size, hopping_distance in reach_cases:
            for delta_option in (0, 1):
                with suite.case(f"{reach_name}, delta={delta_option}"):
                    run_probability_case(
                        executable, grid_size, hopping_distance, delta_option)

        with suite.case("undersized interaction-template reach"):
            run_rejection_case(
                executable, 100, 10, 0.75,
                "Grid size is too small for the requested interaction-distance support.",
                "grasshopper3d-undersized-reach-")

        with suite.case("full grid occupancy"):
            run_rejection_case(
                executable, 1000, 10, 0.1,
                "Number of spins must satisfy 0 < N < grid volume.",
                "grasshopper3d-full-grid-")

        with suite.case("grid volume outside flattened-index range"):
            run_rejection_case(
                executable, 100, 1291, 0.1,
                "Grid size is outside the supported flattened-index range.",
                "grasshopper3d-grid-volume-")

        with suite.case("malformed CLI preserves existing result.dat"):
            run_malformed_cli_preserves_result(executable)

        with suite.case("short positive annealing schedule"):
            run_short_annealing_schedule(executable)

        with suite.case("final partial-round acceptance accounting"):
            run_partial_round_acceptance_case(executable)

        with suite.case("default output rejection"):
            run_default_output_rejection(executable)

        with suite.case("explicit output overwrite"):
            run_explicit_overwrite(executable)

        with suite.case("loaded input preservation during overwrite"):
            run_load_overwrite_preserves_input(executable)

        with suite.case("output cleanup preflight failure"):
            run_output_cleanup_preflight_failure(executable)

        with suite.case("invalid loaded input preserves existing output"):
            run_invalid_load_preserves_existing_output(executable)

        with suite.case("invalid loaded input creates no outputs"):
            run_invalid_load_creates_no_outputs(executable)

        with suite.case("independent final-state energy audit"):
            run_final_state_energy_audit(executable)
    except Exception as error:
        suite.failures.append(("integration test harness", error))

    return suite.finish()


if __name__ == "__main__":
    sys.exit(main())
