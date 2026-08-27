#!/usr/bin/env python3

import contextlib
import math
import pathlib
import subprocess
import sys
import tempfile


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


def require_current_precision(actual, expected, description):
    require(math.isclose(actual, expected, rel_tol=5.0e-6, abs_tol=5.0e-6),
            f"{description}: {actual:.17g} != {expected:.17g}")


def run_probability_case(executable, delta_option):
    total_spins = 100
    grid_size = 10
    hopping_distance = 0.5

    # A centered 4 x 5 x 5 block. For these parameters, the template's
    # four-cell reach contains the full compact support of both kernels.
    coordinates = [
        z * grid_size * grid_size + y * grid_size + x
        for z in range(2, 7)
        for y in range(2, 7)
        for x in range(3, 7)
    ]
    require(len(coordinates) == total_spins,
            "integration fixture has the wrong number of coordinates")

    with tempfile.TemporaryDirectory(
            prefix=f"grasshopper3d-delta{delta_option}-") as directory:
        working_directory = pathlib.Path(directory)
        (working_directory / "initconf.dat").write_text(
            "".join(f"{coordinate}\n" for coordinate in coordinates),
            encoding="utf-8")

        command = [
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
            "-initconf", "load",
            "-delta", str(delta_option),
            "-randomseed", "12345",
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
        reference_probability = reference_energy / (
            2.0 * math.pi * hopping_distance ** 2 * total_spins ** (5.0 / 3.0))

        initial_energy = float(
            (working_directory / "energies.dat")
            .read_text(encoding="utf-8").splitlines()[0])
        reported_energy = float(result_values["final energy"])
        reported_probability = float(result_values["final probability"])

        require_current_precision(
            initial_energy, reference_energy, "initial reported energy")
        require_current_precision(
            reported_energy, reference_energy, "zero-step result energy")
        require_current_precision(
            reported_probability, reference_probability,
            "zero-step result probability")


def main():
    if len(sys.argv) != 2:
        print(f"usage: {pathlib.Path(sys.argv[0]).name} GRASSHOPPER_EXECUTABLE",
              file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    suite = IntegrationTestSuite()

    try:
        require(executable.is_file(), f"executable not found: {executable}")
        for delta_option in (0, 1):
            with suite.case(f"loaded 3D probability with delta={delta_option}"):
                run_probability_case(executable, delta_option)
    except Exception as error:
        suite.failures.append(("integration test harness", error))

    return suite.finish()


if __name__ == "__main__":
    sys.exit(main())
