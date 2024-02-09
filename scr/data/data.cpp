#include "data.h"


std::map<float, float> AirCriterium(std::ifstream& infile) {
	std::map<float, float> air_criterium;
	std::pair<float, float> local_t;
	std::string out;
	std::string temp;
	int i = 0;
	while (std::getline(infile, out)) {
		if (i < 2) {
			i++;
			continue;
		}
		else {
			for (const char& c : out) {
				if ((c != '\t')) {
					temp += c;
				}
				else {
					local_t.first = stod(temp);
					temp.clear();
				}
			}
			local_t.second = stod(temp);
			temp.clear();
			i++;
		}
		air_criterium[local_t.first] = local_t.second;
	}
	return air_criterium;
}
std::map<float, float> Linealize(const std::map<float, float>& unlinealize) {
	std::map<float, float> linealize;
	std::pair<float, float> pair_unlinealize;
	float line_ceil;
	int i = 0;
	for (const auto& [key, value] : unlinealize) {
		if (i == 0) {
			line_ceil = std::ceil(key / 10);
			pair_unlinealize.first = key;
			pair_unlinealize.second = value;
			i++;
		}
		else {
			while (line_ceil < std::ceil(key / 10)) {
				linealize[line_ceil*10] = pair_unlinealize.second + (line_ceil - pair_unlinealize.first / 10) * value;
				line_ceil++;
			}
			pair_unlinealize.first = key;
			pair_unlinealize.second = value;		
		}
	}
	return linealize;
}
//считывает свойство воздуха
void OpenPropertiesOFAir(std::map<float, float>& propertie_of_air, const std::string& file_name) {
	std::ifstream file;
	file.open("data/" + file_name + ".txt", std::ios_base::in);
	if (file.is_open())
		propertie_of_air = Linealize(AirCriterium(file));
	file.close();
}



