#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <array>
#include <string>
#include <vector>
#include <map>

#include <omp.h>
#include "interpolate/interpolate.h"
#include "data/data.h"
#include "read_bmt/bmt_parser.h"
#include "read_begin_conditional.h"
#include "analysis/effective.h"
#include "analysis/Nu.h"
//#include "analysis/Nu.h"

int main() {
	//1,3436*x-9,4642
	std::cin.tie(nullptr);
	std::ios_base::sync_with_stdio(false);
	std::cout << "Enter input file" << std::endl;
	std::string inputFile;
	std::cin >> inputFile;
	std::ifstream infile(inputFile);
	std::string ThermoCharacteristic;
	infile >> ThermoCharacteristic;
	if (ThermoCharacteristic == "Eff") {
		int width_pix = 1280, height_pix = 960, margin_from_edges_pix = 10;
		infile >> width_pix >> height_pix >> margin_from_edges_pix; // 1280 960
		//при использовании стекла, параметры линейной апроксимации температуры
		float linear_ax = 1, linear_b = 0; //1.2436, 9.4642
		//начало координат, обезразмеривание
		float x_begin_coord = 0, y_begin_coord = 0;//0, +29.5 (в мм)
		float non_dimensionalize_coeff = 1; // 4
		infile >> linear_ax >> linear_b >> x_begin_coord >> y_begin_coord >> non_dimensionalize_coeff;
		std::string x_name = "z/d", y_name = "x/d", temp_name = "T", eff_name = "eff";
		infile >> x_name >> y_name >> temp_name >> eff_name;
		std::string beginParametr = "table_in_c.txt"; //файл с начальными параметрами
		infile >> beginParametr;
		std::vector<beginConditionals> bC = InputBeginConditionals(beginParametr);
		std::vector<std::string> parametr_calculate;
		
		while (!infile.eof()) {
			
			std::string parcal;
			infile >> parcal;
			parametr_calculate.push_back(parcal);
		}
		infile.close();
		int i = 0;
		while(i < parametr_calculate.size()) {
			
			if (parametr_calculate[i] == "AREA_LOCAL") {
				AreaEffectiveness(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, temp_name, eff_name);
				i++;
			}
			else if (parametr_calculate[i] == "AREA_INTEGRAL") {
				AreaIntegralEffectiveness(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, temp_name, eff_name);
				i++;
			}
			else if (parametr_calculate[i] == "ALONG_X_LOCAL_xb") {
				AlongXbEffectiveness(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, temp_name, eff_name);
				i++;		
			}
			else if (parametr_calculate[i] == "ALONG_X_LOCAL") {
				AlongXEffectiveness(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, temp_name, eff_name, std::stof(parametr_calculate[i+1]));
				i += 2;
			}
			else if (parametr_calculate[i] == "ALONG_Y_LOCAL") {
				AlongYEffectiveness(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, temp_name, eff_name, std::stof(parametr_calculate[i + 1]));
				i += 2;

			}
			else if (parametr_calculate[i] == "ALONG_X_INTEGRAL_xb") {
				AlongXbIntegralEffectiveness(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, temp_name, eff_name, std::stof(parametr_calculate[i + 1]));
				i += 2;
			}
			else if (parametr_calculate[i] == "ALONG_X_INTEGRAL") {
				AlongXIntegralEffectiveness(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, temp_name, eff_name, std::stof(parametr_calculate[i + 1]), std::stof(parametr_calculate[i + 2]));
				i += 3;

			}
			else if (parametr_calculate[i] == "ALONG_Y_INTEGRAL") {
				AlongYIntegralEffectiveness(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, temp_name, eff_name, std::stof(parametr_calculate[i + 1]), std::stof(parametr_calculate[i + 2]));
				i += 3;
			}
			else {
				std::cout << "uncorrect parametr read" << std::endl;
				break;
			}
		}
	}
	else if (ThermoCharacteristic == "Nu4") {
		
		std::cout << "Nu4\n";
		int width_pix = 1280, height_pix = 960, margin_from_edges_pix = 10;
		infile >> width_pix >> height_pix >> margin_from_edges_pix; // 1280 960
		//при использовании стекла, параметры линейной апроксимации температуры
		float linear_ax = 1, linear_b = 0; //1.2436, 9.4642
		//начало координат, обезразмеривание
		float x_begin_coord = 0, y_begin_coord = 0;//0, +29.5 (в мм)
		float non_dimensionalize_coeff = 1; // 4
		infile >> linear_ax >> linear_b >> x_begin_coord >> y_begin_coord >> non_dimensionalize_coeff;
		float lenght_plate = 0.520;
		float weight_plate = 0.382;
		infile >> lenght_plate >> weight_plate;
		float epsilon_kraska, epsilon_material;
		infile >> epsilon_kraska  >> epsilon_material;
		float diameter_jet;
		infile >> diameter_jet;
		std::string x_name = "z/d", y_name = "x/d";
		infile >> x_name >> y_name;
		std::string beginParametr = "table_in_c.txt"; //файл с начальными параметрами
		infile >> beginParametr;
		std::vector<beginConditionalsNu4> bC = InputBeginConditionalsNu4(beginParametr);
		std::vector<std::string> parametr_calculate;
		while (!infile.eof()) {

			std::string parcal;
			infile >> parcal;
			parametr_calculate.push_back(parcal);
		}
		infile.close();
		for (int j = 0; j < bC.size(); j++) {
			int i = 0;
			Nu object(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, j, lenght_plate, weight_plate, epsilon_kraska, epsilon_material, diameter_jet);
			while (i < parametr_calculate.size()) {

				if (parametr_calculate[i] == "AREA_LOCAL") {
					object.AreaNu();
					i++;
					std::cout << bC[j].number << std::endl;
				}
				else if (parametr_calculate[i] == "AREA_INTEGRAL") {
					object.AreaIntegralNu();
					i++;
				}
				else if (parametr_calculate[i] == "ALONG_X_LOCAL_xb") {
					object.AlongXbNu();
					i++;
				}
				else if (parametr_calculate[i] == "ALONG_X_LOCAL") {
					object.AlongXNu(std::stof(parametr_calculate[i + 1]));
					i += 2;
				}
				else if (parametr_calculate[i] == "ALONG_Y_LOCAL") {
					object.AlongYNu(std::stof(parametr_calculate[i + 1]));
					i += 2;

				}
				else if (parametr_calculate[i] == "ALONG_X_INTEGRAL_xb") {
					object.AlongXbIntegralNu(std::stof(parametr_calculate[i + 1]));
					i += 2;
				}
				else if (parametr_calculate[i] == "ALONG_X_INTEGRAL") {
					object.AlongXIntegralNu(std::stof(parametr_calculate[i + 1]), std::stof(parametr_calculate[i + 2]));
					i += 3;

				}
				else if (parametr_calculate[i] == "ALONG_Y_INTEGRAL") {
					object.AlongYIntegralNu(std::stof(parametr_calculate[i + 1]), std::stof(parametr_calculate[i + 2]));
					i += 3;
				}
				else if (parametr_calculate[i] == "ROUND_INTEGRAL_Nu") {
					object.RoundIntegralNu(std::stof(parametr_calculate[i + 1]), std::stof(parametr_calculate[i + 2]), std::stof(parametr_calculate[i + 3]));
					i += 4;
				}
				else {
					std::cout << "uncorrect parametr read" << std::endl;
					break;
				}
			}
		}
	}
	else if (ThermoCharacteristic == "Nu3") {
		std::cout << "Nu3\n";
		int width_pix = 1280, height_pix = 960, margin_from_edges_pix = 10;
		infile >> width_pix >> height_pix >> margin_from_edges_pix; // 1280 960
		//при использовании стекла, параметры линейной апроксимации температуры
		float linear_ax = 1, linear_b = 0; //1.2436, 9.4642
		//начало координат, обезразмеривание
		float x_begin_coord = 0, y_begin_coord = 0;//0, +29.5 (в мм)
		float non_dimensionalize_coeff = 1; // 4
		infile >> linear_ax >> linear_b >> x_begin_coord >> y_begin_coord >> non_dimensionalize_coeff;
		float lenght_plate = 0.520;
		float weight_plate = 0.382;
		infile >> lenght_plate >> weight_plate;
		float epsilon_kraska, epsilon_material;
		infile >> epsilon_kraska >> epsilon_material;
		float diameter_jet;
		infile >> diameter_jet;
		std::string x_name = "z/d", y_name = "x/d";
		infile >> x_name >> y_name;
		std::string beginParametr = "table_in_c.txt"; //файл с начальными параметрами
		infile >> beginParametr;
		std::vector<beginConditionalsNu4> bC = InputBeginConditionalsNu4(beginParametr);
		std::vector<std::string> parametr_calculate;

		while (!infile.eof()) {

			std::string parcal;
			infile >> parcal;
			parametr_calculate.push_back(parcal);
		}
		infile.close();
		int i = 0;
		for (int j = 0; j < bC.size(); j++) {
			Nu object(width_pix, height_pix, margin_from_edges_pix, linear_ax, linear_b, x_begin_coord, y_begin_coord, non_dimensionalize_coeff, bC, x_name, y_name, j, lenght_plate, weight_plate, epsilon_kraska, epsilon_material, diameter_jet);
			while (i < parametr_calculate.size()) {

				if (parametr_calculate[i] == "AREA_LOCAL") {
					object.AreaNu();
					i++;
				}
				else if (parametr_calculate[i] == "AREA_INTEGRAL") {
					object.AreaIntegralNu();
					i++;
				}
				else if (parametr_calculate[i] == "ALONG_X_LOCAL_xb") {
					object.AlongXbNu();
					i++;
				}
				else if (parametr_calculate[i] == "ALONG_X_LOCAL") {
					object.AlongXNu(std::stof(parametr_calculate[i + 1]));
					i += 2;
				}
				else if (parametr_calculate[i] == "ALONG_Y_LOCAL") {
					object.AlongYNu(std::stof(parametr_calculate[i + 1]));
					i += 2;

				}
				else if (parametr_calculate[i] == "ALONG_X_INTEGRAL_xb") {
					object.AlongXbIntegralNu(std::stof(parametr_calculate[i + 1]));
					i += 2;
				}
				else if (parametr_calculate[i] == "ALONG_X_INTEGRAL") {
					object.AlongXIntegralNu(std::stof(parametr_calculate[i + 1]), std::stof(parametr_calculate[i + 2]));
					i += 3;

				}
				else if (parametr_calculate[i] == "ALONG_Y_INTEGRAL") {
					object.AlongYIntegralNu(std::stof(parametr_calculate[i + 1]), std::stof(parametr_calculate[i + 2]));
					i += 3;
				}
				else if (parametr_calculate[i] == "ROUND_INTEGRAL_Nu") {
					object.RoundIntegralNu(std::stof(parametr_calculate[i + 1]), std::stof(parametr_calculate[i + 2]), std::stof(parametr_calculate[i + 3]));
					i += 4;
				}
				else {
					std::cout << "uncorrect parametr read" << std::endl;
					break;
				}
			}
		}
	}
	else if (ThermoCharacteristic == "BMT") {

		std::string beginParametr = "1.BMT"; //файл с начальными параметрами
		infile >> beginParametr;
		DataExtract a(beginParametr);
		a.DataPrintGraphTable();
		a.DataPrintImage();
		//std::vector<float> x = a.DataExportTemperature();
		//std::vector<Color> y = a.DataExportImage();
	}
	else {
		std::cout << "Error\n";
	}
	return 0;
}

//std::cin.tie(nullptr);
		//std::ios_base::sync_with_stdio(false);
		//float T = 34.5;
		//std::map<float, float> AirKinematicViscosity;
		//PROPERTIES_OF_AIR(AirKinematicViscosity);
		//std::cout << LinearInterpolation(AirKinematicViscosity, T) << std::endl;