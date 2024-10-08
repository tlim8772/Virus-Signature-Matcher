#pragma once

struct MatchResult {
    std::string sample_name;
    std::string signature_name;
    double match_score;
};

void runMatcher(const std::vector<klibpp::KSeq>& samples,
                const std::vector<klibpp::KSeq>& signatures,
                std::vector<MatchResult>& matches);
