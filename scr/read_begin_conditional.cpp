#include "read_begin_conditional.h"
std::vector<beginConditionals> InputBeginConditionals(const std::string& nameBeginConditionalFile) {
	std::ifstream infile(nameBeginConditionalFile);
	// Temporary buffer
	std::string out;
	std::vector<beginConditionals> bC;
	while (std::getline(infile, out)) {
		int i = 0;
		std::vector<std::string> inp(18);
		for (const char& c : out) {
			if (c != '\t') {
				inp[i] += c;
			}
			else {
				i++;
			}
		}
		bC.push_back({ "IR00" + inp[0],
			std::stof(inp[1]),
			std::stof(inp[2]),
			std::stof(inp[3]),
			std::stof(inp[4]),
			std::stof(inp[5]),
			std::stof(inp[6]),
			std::stof(inp[7]),
			std::stoi(inp[8]) * 2,
			std::stoi(inp[9]) * 2,
			std::stoi(inp[10]) * 2,
			std::stoi(inp[11]) * 2,
			std::stoi(inp[12]) * 2,
			std::stoi(inp[13]) * 2,
			std::stoi(inp[14]) * 2,
			std::stoi(inp[15]) * 2,
			std::stof(inp[16]) / 2,
			});
	}
	return bC;
}