// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Olga Goulko and David Llamas

#include "doctest/doctest.h"
#include "output.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <system_error>

namespace {

class ScopedTemporaryPath {
public:
    explicit ScopedTemporaryPath(const std::string& label)
        : path_(std::filesystem::temp_directory_path()
                /("grasshopper3d-"+label+"-"
                  +std::to_string(
                      std::chrono::steady_clock::now().time_since_epoch().count()))) {}

    ~ScopedTemporaryPath()
        {
        std::error_code error;
        std::filesystem::remove_all(path_,error);
        }

    const std::filesystem::path& path() const {return path_;}

private:
    std::filesystem::path path_;
};

template<typename Operation>
void checkRuntimeErrorContains(Operation operation, const std::string& expectedText)
    {
    try
        {
        operation();
        FAIL_CHECK("output failure was not reported");
        }
    catch(const std::runtime_error& error)
        {
        CHECK(std::string(error.what()).find(expectedText)!=std::string::npos);
        }
    }

bool devFullAvailable()
    {
#ifdef __linux__
    std::ofstream probe("/dev/full");
    return probe.is_open();
#else
    return false;
#endif
    }

} // namespace

TEST_CASE("checked output opens writes flushes and closes a file")
    {
    ScopedTemporaryPath temporaryPath("checked-output");
    const std::string filename=temporaryPath.path().string();
    std::ofstream output;
    openOutputFile(output,filename);
    output << "checked output\n";
    finishOutputFile(output,filename);

    CHECK_FALSE(output.is_open());
    std::ifstream input(filename);
    REQUIRE(input.is_open());
    const std::string contents((std::istreambuf_iterator<char>(input)),
                               std::istreambuf_iterator<char>());
    CHECK(contents=="checked output\n");
    }

TEST_CASE("bounded configuration output stops snapshots without failing the run")
    {
    ScopedTemporaryPath temporaryPath("bounded-output");
    const std::string filename=temporaryPath.path().string();
    BoundedOutputFile output(filename,10);

    CHECK(output.writeIfFits("12345\n"));
    CHECK_FALSE(output.writeIfFits("abcdef\n"));
    CHECK(output.limitReached());
    CHECK_FALSE(output.writeIfFits("x\n"));
    CHECK(output.bytesWritten()==6);
    output.finish();

    CHECK(std::filesystem::file_size(filename)<=10);
    std::ifstream input(filename);
    REQUIRE(input.is_open());
    const std::string contents((std::istreambuf_iterator<char>(input)),
                               std::istreambuf_iterator<char>());
    CHECK(contents=="12345\n");
    }

TEST_CASE("output open failures identify the target file")
    {
    ScopedTemporaryPath temporaryPath("missing-output-parent");
    const std::string filename=(temporaryPath.path()/"missing"/"output.dat").string();
    checkRuntimeErrorContains([&filename]()
        {
        std::ofstream output;
        openOutputFile(output,filename);
        },filename);
    }

TEST_CASE("post-open write failures identify dev full"
          *doctest::skip(!devFullAvailable()))
    {
    std::ofstream output;
    openOutputFile(output,"/dev/full");
    output << "data\n";
    checkRuntimeErrorContains([&output]()
        {
        finishOutputFile(output,"/dev/full");
        },"/dev/full");
    }
