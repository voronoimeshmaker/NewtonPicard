#pragma once

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>


namespace study {

class TeeBuffer final : public std::streambuf {
public:
    TeeBuffer() = default;

    void set_targets(std::streambuf* first, std::streambuf* second) noexcept {
        first_ = first;
        second_ = second;
    }

protected:
    int overflow(int ch) override {
        if (traits_type::eq_int_type(ch, traits_type::eof())) {
            return traits_type::not_eof(ch);
        }

        const char c = traits_type::to_char_type(ch);
        const auto r1 = first_ != nullptr ? first_->sputc(c) : ch;
        const auto r2 = second_ != nullptr ? second_->sputc(c) : ch;

        if ((first_ != nullptr && traits_type::eq_int_type(r1, traits_type::eof())) ||
            (second_ != nullptr && traits_type::eq_int_type(r2, traits_type::eof()))) {
            return traits_type::eof();
        }

        return ch;
    }

    std::streamsize xsputn(const char* s, std::streamsize count) override {
        const auto n1 = first_ != nullptr ? first_->sputn(s, count) : count;
        const auto n2 = second_ != nullptr ? second_->sputn(s, count) : count;
        return (n1 == count && n2 == count) ? count : 0;
    }

    int sync() override {
        const int r1 = first_ != nullptr ? first_->pubsync() : 0;
        const int r2 = second_ != nullptr ? second_->pubsync() : 0;
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }

private:
    std::streambuf* first_{nullptr};
    std::streambuf* second_{nullptr};
};

/// CHANGED: the default console log no longer truncates.
///
/// Thirteen studies construct this class with no argument, which opened
/// Result/console.log with std::ios::trunc.  Every run of any of them erased
/// the log of the previous one, and of every other study in that directory, so
/// the protocol of keeping all runs could not be honoured through this path at
/// all.  The two constructors now say which case is which:
///
///   ConsoleTee()          shared Result/console.log, opened for append and
///                         preceded by a separator naming the UTC start time.
///                         Nothing is ever lost.
///   ConsoleTee(path)      a log inside a directory the caller owns, opened
///                         for truncation.  A study that calls
///                         create_result_directory() first gets a fresh
///                         directory per run, so truncation there is correct
///                         and this is the form to prefer.
class ConsoleTee final {
public:
    /// Shared log, appended to.
    ConsoleTee() : ConsoleTee(default_path(), false) {}

    /// A log the caller names.  Truncates, because the caller is expected to
    /// own the directory; pass a per-run directory rather than a fixed name.
    explicit ConsoleTee(const std::filesystem::path& path) : ConsoleTee(path, true) {}

    ConsoleTee(const std::filesystem::path& path, bool truncate)
    {
        const auto parent = path.parent_path();
        if (!parent.empty()) {
            std::filesystem::create_directories(parent);
        }

        log_.open(path, truncate ? (std::ios::out | std::ios::trunc)
                                 : (std::ios::out | std::ios::app));
        if (!log_) {
            throw std::runtime_error(
                "Could not open study console log: " + path.string());
        }

        if (!truncate) {
            log_ << "\n===== run started " << utc_stamp() << " =====\n";
            log_.flush();
        }

        old_cout_ = std::cout.rdbuf();
        old_cerr_ = std::cerr.rdbuf();

        cout_tee_.set_targets(old_cout_, log_.rdbuf());
        cerr_tee_.set_targets(old_cerr_, log_.rdbuf());

        std::cout.rdbuf(&cout_tee_);
        std::cerr.rdbuf(&cerr_tee_);
    }

    ConsoleTee(const ConsoleTee&) = delete;
    ConsoleTee& operator=(const ConsoleTee&) = delete;

    ~ConsoleTee() {
        std::cout.flush();
        std::cerr.flush();

        if (old_cout_ != nullptr) {
            std::cout.rdbuf(old_cout_);
        }
        if (old_cerr_ != nullptr) {
            std::cerr.rdbuf(old_cerr_);
        }
    }

private:
    [[nodiscard]] static std::filesystem::path default_path() {
        return std::filesystem::path{"Result"} / "console.log";
    }

    [[nodiscard]] static std::string utc_stamp() {
        const auto now = std::chrono::system_clock::now();
        const std::time_t raw = std::chrono::system_clock::to_time_t(now);
        std::tm value{};
#if defined(_WIN32)
        gmtime_s(&value, &raw);
#else
        gmtime_r(&raw, &value);
#endif
        std::ostringstream out;
        out << std::put_time(&value, "%Y-%m-%dT%H:%M:%SZ");
        return out.str();
    }

    std::ofstream log_;
    std::streambuf* old_cout_{nullptr};
    std::streambuf* old_cerr_{nullptr};
    TeeBuffer cout_tee_;
    TeeBuffer cerr_tee_;
};

} // namespace study
