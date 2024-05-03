#include "Nu.h"


Nu::Nu(std::string Files, int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, std::vector<beginConditionalsNu4>& bC, const std::string& x_name, const std::string& y_name, int i, float lenght_plate, float weight_plate, float epsilon_kraska, float epsilon_material, float hydravlic_diameter, const std::string& folder)
: width_pix(width_pix), height_pix(height_pix), margin_from_edges_pix(margin_from_edges_pix), linear_ax(linear_ax), linear_b(linear_b), y_begin_coord(y_begin_coord), x_begin_coord(x_begin_coord),non_dimensionalize_coeff(non_dimensionalize_coeff), bC(bC), x_name(x_name), y_name(y_name), i(i), folder(folder)
{
	//number expiriment	numjets	s/d	h/d	Q,g/s	U	A	W	tair	tjet	x0	y0	ymin	ymax dx 
	//bC[i].dx_mm
	
	PROPERTIES_OF_AIR(AirDensity);
	PROPERTIES_OF_AIR(AirSpecificHeat);
	
	PROPERTIES_OF_AIR(AirPrandtl);
	PROPERTIES_OF_AIR(AirThermalDiffusivity);
	PROPERTIES_OF_AIR(AirKinematicViscosity);
	PROPERTIES_OF_AIR(AirLambdaPatametr);
	PROPERTIES_OF_AIR(AirVolumetric—oefficientofThermalExpansion);
	int k = 0;
	float Ti;
	float AirLambdaTair = LinearInterpolation(AirLambdaPatametr, bC[i].Tair);
	float AirPrandtlTair = LinearInterpolation(AirPrandtl, bC[i].Tair);
	float AirKinematicViscosityTair = LinearInterpolation(AirKinematicViscosity, bC[i].Tair);
	float AirVolumetric—oefficientofThermalExpansionTair = LinearInterpolation(AirVolumetric—oefficientofThermalExpansion, bC[i].Tair);
	float AirThermalDiffusivityTair = LinearInterpolation(AirThermalDiffusivity, bC[i].Tair);
    std::ifstream infile(Files + "/" + bC[i].number + ".txt");
    for (int j = 0; j < height_pix * width_pix; j++) {
        infile >> Ti;
		if ((j % width_pix + margin_from_edges_pix < bC[i].xmax) && (j % width_pix - margin_from_edges_pix > bC[i].xmin) && (j / width_pix + margin_from_edges_pix < bC[i].ymax) && (j / width_pix - margin_from_edges_pix > bC[i].ymin)) {
#pragma omp parallel
			{
				T_.push_back(linear_ax * Ti - linear_b );//ÔÓ„Â¯ÌÓÒÚ¸ ÒÚÂÍÎ‡
				 //
				float AirLambdaTk = LinearInterpolation(AirLambdaPatametr, T_[k]);
				float AirPrandtlTk = LinearInterpolation(AirPrandtl, T_[k]);
				X.push_back(((j % width_pix) * bC[i].dx_mm - bC[i].x0 * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
				Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
				Qizl_.push_back(sigma * (epsilon_kraska + epsilon_material) * (pow((T_[k] + 273.15), 4) - pow((273.15 + bC[i].Tair), 4)));
				Qel_.push_back(bC[i].W / (weight_plate * lenght_plate));
				//Ra_.push_back(LinearInterpolation(AirVolumetric—oefficientofThermalExpansion, T_[k]) * g_free * (Qel_[k] - Qizl_[k]) * pow(diameter_jet, 4) / AirLambdaTk / LinearInterpolation(AirKinematicViscosity, T_[k]) / LinearInterpolation(AirThermalDiffusivity, T_[k]));
				Ra2_.push_back(AirVolumetric—oefficientofThermalExpansionTair * g_free * (Qel_[k] - Qizl_[k]) * pow(hydravlic_diameter, 4) / AirLambdaTair / AirKinematicViscosityTair / AirThermalDiffusivityTair);
				//Qck1_.push_back(0.508 * AirLambdaTk * (T_[k] - bC[i].Tair) * pow((Ra_[k] * AirPrandtlTk / (0.952 + AirPrandtlTk)), 0.25) / diameter_jet);
				Qck2_.push_back(0.508 * AirLambdaTair * (T_[k] - bC[i].Tair) * pow((Ra2_[k] * AirPrandtlTair / (0.952 + AirPrandtlTair)), 0.25) / hydravlic_diameter);
				Qw_.push_back((Qel_[k] - Qizl_[k] - Qck2_[k]) >= 0 ? (Qel_[k] - Qizl_[k] - Qck2_[k]) : 0);
				Nu_.push_back(hydravlic_diameter * (Qw_[k]) / (T_[k] - bC[i].Tjet ) / AirLambdaTk); //- 1.2 * (80 - bC[i].Tjet) / 80
				Alpha_.push_back(Qw_[k] / (T_[k] - bC[i].Tjet));
				k++;
			}
        }
    }
    infile.close();
}

Nu::Nu(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, std::vector<beginConditionalsNu3>& bC, const std::string& x_name, const std::string& y_name, int i, float lenght_plate, float weight_plate, float epsilon_kraska, float epsilon_material, float hydravlic_diameter, int number_of_wall, float x_begin, float high, float lenght, float weight, const std::string& folder)
{
	PROPERTIES_OF_AIR(AirDensity);
	PROPERTIES_OF_AIR(AirSpecificHeat);
	PROPERTIES_OF_AIR(AirPrandtl);
	PROPERTIES_OF_AIR(AirKinematicViscosity);
	PROPERTIES_OF_AIR(AirLambdaPatametr);
	int k = 0;
	float Ti;
	std::ifstream infile(bC[i].number + ".txt");
	for (int j = 0; j < height_pix * width_pix; j++) {
		infile >> Ti;
		if ((j % width_pix + margin_from_edges_pix < bC[i].xmax) && (j % width_pix - margin_from_edges_pix > bC[i].xmin) && (j / width_pix - margin_from_edges_pix < bC[i].ymax) && (j / width_pix + margin_from_edges_pix > bC[i].ymin)) {
#pragma omp parallel
			{
				T_.push_back(linear_ax * Ti - linear_b);
				X.push_back(((j % width_pix) * bC[i].dx_mm - bC[i].x0 * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
				Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
				Tin_.push_back((bC[i].Tair + bC[i].W * (x_begin - X[k])) / (LinearInterpolation(AirSpecificHeat, T_[k]) * LinearInterpolation(AirDensity, T_[k]) * bC[i].velosity * high * lenght * weight));
				float Tin0 = T_[k] - 1 / 3 * Tin_[k];
				float Tin1 = T_[k] - 1 / 3 * bC[i].Tair;
				float AirPrandtlT0 = LinearInterpolation(AirPrandtl, Tin0);
				float AirPrandtlT1 = LinearInterpolation(AirPrandtl, Tin1);
				
				Qizl_.push_back(sigma * epsilon_kraska  * (pow((T_[k] + 273.15), 4) - pow((273.15 + bC[i].Tair), 4)));
				Qel_.push_back(bC[i].W / (number_of_wall * weight_plate * lenght_plate));
				Ra_.push_back(2.0 * 9.8 * (T_[k] - Tin_[k]) * pow((X[k]), 3) / (T_[k] + Tin_[k] + 2.0 * 273.15) / LinearInterpolation(AirKinematicViscosity, Tin0) / LinearInterpolation(AirThermalDiffusivity, Tin0));
				Ra2_.push_back(2.0 * 9.8 * (T_[k] - bC[i].Tair) * pow((X[k]), 3) / (T_[k] + bC[i].Tair + 2.0 * 273.15) / LinearInterpolation(AirKinematicViscosity, Tin1) / LinearInterpolation(AirThermalDiffusivity, Tin1));
				Qck1_.push_back(0.508 * LinearInterpolation(AirLambdaPatametr, Tin0) * (T_[k] - Tin_[k]) * pow((Ra_[k] * AirPrandtlT0 / (0.952 + AirPrandtlT0)), 0.25) / X[k]);
				Qck2_.push_back(0.508 * LinearInterpolation(AirLambdaPatametr, Tin1) * (T_[k] - bC[i].Tair) * pow((Ra2_[k] * AirPrandtlT1 / (0.952 + AirPrandtlT1)), 0.25) / X[k]);
				Qw_.push_back((Qel_[k] - Qizl_[k] - Qck2_[k]) >= 0 ? (Qel_[k] - Qizl_[k] - Qck2_[k]) : 0);
				Nu_.push_back(hydravlic_diameter * (Qw_[k]) / (T_[k] - bC[i].Tjet) / AirPrandtlT0);
				Alpha_.push_back(Qw_[k] / (T_[k] - bC[i].Tjet));
				k++;
			}
		}
	}
	infile.close();
}

void Nu::AreaNu()
{
	std::string _filename = folder + "/Nu_area_" + bC[i].number + ".txt";

	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << x_name << ", " << y_name << ", " << "T" << ", " << "Qizl" << ", " << "Qel" << ", " << "Ra_out" << ", " << "Qck_out" << ", " << "Alpha" << ", " << "Nu" << '\n';
#pragma omp parallel for
	for (int j = 0; j < T_.size(); j++) {
		f << X[j] << ", " << Y[j] << ", " << T_[j] << ", " << Qizl_[j] << ", " << Qel_[j] <<  ", " << Ra2_[j] << ", "  << Qck2_[j] << ", " << Alpha_[j] << ", " << Nu_[j] << '\n';
	}
	f.close();
}

void Nu::TminNumax(std::map<int, std::pair<float, float>>& a)
{
	//int x0;
	//int y0;
	float Nu_max = 0;
	float T_min = FLT_MAX;
	for (int j = 0; j < T_.size(); j++) {
		if (Nu_[j] > Nu_max) {
			Nu_max = Nu_[j];
			///x0 = j % width_pix;
			//y0 = j / width_pix;
		}
		if (T_[j] < T_min) {
			T_min = T_[j];
		}
	}
	a[i].first = T_min;
	a[i].second = Nu_max;
	//bC[i].x0 = x0;
	//bC[i].y0 = y0;
	std::cout << "i = " << i << "T_min = " << T_min << "Nu_max = " << Nu_max << std::endl;
}

void Nu::AreaIntegralNu()
{
	std::string _filename = folder + "/Nu_area_integral" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	float Nu_integral = 0;
	for (int j = 0; j < T_.size(); j++) {
		Nu_integral += Nu_[j];
	}
	f << "Nu = " << Nu_integral / T_.size();	
	f.close();
}

void Nu::AlongXbNu()
{
	std::string _filename = folder + "/Nu_alongXb_" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << x_name << "\t" << "T" << "\t" << "Qizl" << "\t" << "Qel"  << "\t" << "Ra_out" << "\t" << "Qck_out" << "\t" << "Alpha" << "\t" << "Nu" << '\n';
#pragma omp parallel for
	for (int j = 0; j < T_.size(); j++) {
		if ((Y[j] == 0)) {
			f << X[j] << "\t" << T_[j] << "\t" << Qizl_[j] << "\t" << Qel_[j]  << "\t" << Ra2_[j] << "\t" << Qck2_[j] << "\t" << Alpha_[j] << "\t" << Nu_[j] << '\n';
		}
	}
	f.close();
}

void Nu::AlongXNu(float x_0)
{
	std::string _filename = folder + "/Nu_alongX_" + std::to_string(static_cast<int>(x_0)) + "mm" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << x_name << "\t" << "T" << "\t" << "Qizl" << "\t" << "Qel" << "\t" << "Ra_out" << "\t" << "Qck_out" << "\t" << "Alpha" << "\t" << "Nu" << '\n';
#pragma omp parallel for
	for (int j = 0; j < T_.size(); j++) {
		if (static_cast<int>(std::abs((Y[j] * non_dimensionalize_coeff + x_0) / bC[i].dx_mm) + 0.1) == 0) {
			f << X[j] << "\t" << T_[j] << "\t" << Qizl_[j] << "\t" << Qel_[j] << "\t" << Ra2_[j] << "\t" << Qck2_[j] << "\t" << Alpha_[j] << "\t" << Nu_[j] << '\n';
		}
	}
	f.close();
}

void Nu::AlongYNu(float y_0)
{
	std::string _filename = folder + "/Nu_alongY_" + std::to_string(static_cast<int>(y_0)) + "mm" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << y_name << "\t" << "T" << "\t" << "Qizl" << "\t" << "Qel" << "\t" << "Ra_out" << "\t" << "Qck_out" << "\t" << "Alpha" << "\t" << "Nu" << '\n';
#pragma omp parallel for
	for (int j = 0; j < T_.size(); j++) {
		int we = static_cast<int>((X[j] * non_dimensionalize_coeff + y_0) / bC[i].dx_mm);
		if (static_cast<int>(std::abs((X[j] * non_dimensionalize_coeff + y_0) / bC[i].dx_mm) + 0.1) == 0) {
			f << Y[j] << "\t" << T_[j] << "\t" << Qizl_[j] << "\t" << Qel_[j] << "\t" << Ra2_[j] << "\t" << Qck2_[j] << "\t" << Alpha_[j] << "\t" << Nu_[j] << '\n';
		}
	}
	f.close();
}

void Nu::AlongXbIntegralNu(float dx)
{
	float a = std::min(std::min(std::abs(bC[i].xmax - bC[i].x0), std::abs(bC[i].xmin - bC[i].x0)), std::min(std::abs(bC[i].ymax - bC[i].y0), std::abs(bC[i].ymin - bC[i].y0))) / dx;
	std::string _filename = folder + "/Nu_integtal_Xb_" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << y_name << "\t" << "T" << "\t" << "Qizl" << "\t" << "Qel" << "\t" << "Ra_out" << "\t" << "Qck_out" << "\t" << "Alpha" << "\t" << "Nu" << '\n';
#pragma omp parallel for
	for (float j = 1; j < a; j++) {
		float k = 0;
		float Tint = 0;
		float Qizlint = 0;
		float Qelint = 0;
		float Ra_outint = 0;
		float Qck_outint = 0;
		float Alpha_int = 0;
		float Nu_int = 0;
		for (int t = 0; t < T_.size(); t++) {
			if ((std::abs(X[t]) <= j * dx / non_dimensionalize_coeff) && (Y[t] == 0)) {
				Tint += T_[t];
				Qizlint += Qizl_[t];
				Qelint += Qel_[t];
				Ra_outint += Ra2_[t];
				Qck_outint += Qck2_[t];
				Alpha_int += Alpha_[t];
				Nu_int += Nu_[t];
				k++;
			}
		}
		if (k > 0) {
			f << j * dx / non_dimensionalize_coeff << "\t" << Tint / k << "\t" << Qizlint / k << "\t" << Qelint / k << "\t" << Ra_outint / k << "\t" << Qck_outint / k << "\t" << Alpha_int / k << "\t" << Nu_int / k << '\n';
		}
	}
	f.close();
}

void Nu::AlongXIntegralNu(float dx, float x_0)
{
	float a = std::min(std::min(std::abs(bC[i].xmax - bC[i].x0 - x_0 / non_dimensionalize_coeff), std::abs(bC[i].xmin - bC[i].x0 - x_0 / non_dimensionalize_coeff)), static_cast<float>(std::min(std::abs(bC[i].ymax - bC[i].y0), std::abs(bC[i].ymin - bC[i].y0)))) / dx;
	std::string _filename = folder + "/Nu_integtal_x0_" + std::to_string(static_cast<int>(x_0)) + "_" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << y_name << "\t" << "T" << "\t" << "Qizl" << "\t" << "Qel" << "\t" << "Ra_out" << "\t" << "Qck_out" << "\t" << "Alpha" << "\t" << "Nu" << '\n';
#pragma omp parallel for
	for (float j = 1; j < a; j++) {
		float k = 0;
		float Tint = 0;
		float Qizlint = 0;
		float Qelint = 0;
		float Ra_outint = 0;
		float Qck_outint = 0;
		float Alpha_int = 0;
		float Nu_int = 0;
		for (int t = 0; t < T_.size(); t++) {
			if ((std::abs(X[t]) <= j * dx / non_dimensionalize_coeff) && (static_cast<int>((Y[t] * non_dimensionalize_coeff + x_0) / bC[i].dx_mm) == 0)){
				Tint += T_[t];
				Qizlint += Qizl_[t];
				Qelint += Qel_[t];
				Ra_outint += Ra2_[t];
				Qck_outint += Qck2_[t];
				Alpha_int += Alpha_[t];
				Nu_int += Nu_[t];
				k++;
			}
		}
		if (k > 0) {
			f << j * dx / non_dimensionalize_coeff << "\t" << Tint / k << "\t" << Qizlint / k << "\t" << Qelint / k << "\t" << Ra_outint / k << "\t" << Qck_outint / k << "\t" << Alpha_int / k << "\t" << Nu_int / k << '\n';
		}
	}
	f.close();
}

void Nu::AlongYIntegralNu(float dx, float y_0)
{
	float a = std::min(static_cast<float>(std::min(std::abs(bC[i].xmax - bC[i].x0), std::abs(bC[i].xmin - bC[i].x0))), std::min(std::abs(bC[i].ymax - bC[i].y0 - y_0 / non_dimensionalize_coeff), std::abs(bC[i].ymin - bC[i].y0 - y_0 / non_dimensionalize_coeff))) / dx;
	std::string _filename = folder + "/Nu_integtal_y0_" + std::to_string(static_cast<int>(y_0)) + "_" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << y_name << "\t" << "T" << "\t" << "Qizl" << "\t" << "Qel" << "\t" << "Ra_out" << "\t" << "Qck_out" << "\t" << "Alpha" << "\t" << "Nu" << '\n';
#pragma omp parallel for
	for (float j = 1; j < a; j++) {
		float k = 0;
		float Tint = 0;
		float Qizlint = 0;
		float Qelint = 0;
		float Ra_outint = 0;
		float Qck_outint = 0;
		float Alpha_int = 0;
		float Nu_int = 0;
		for (int t = 0; t < T_.size(); t++) {
			if ((std::abs(Y[t]) <= j * dx / non_dimensionalize_coeff) && (static_cast<int>((X[t] * non_dimensionalize_coeff + y_0) / bC[i].dx_mm) == 0)) {
				Tint += T_[t];
				Qizlint += Qizl_[t];
				Qelint += Qel_[t];
				Ra_outint += Ra2_[t];
				Qck_outint += Qck2_[t];
				Alpha_int += Alpha_[t];
				Nu_int += Nu_[t];
				k++;
			}
		}
		if (k > 0) {
			f << j * dx / non_dimensionalize_coeff << "\t" << Tint / k << "\t" << Qizlint / k << "\t" << Qelint / k << "\t" << Ra_outint / k << "\t" << Qck_outint / k << "\t" << Alpha_int / k << "\t" << Nu_int / k << '\n';
		}
	}
	f.close();
}

void Nu::RoundIntegralNu(float x_0, float y_0, float dx)
{
	float a = std::min(std::min(std::abs(bC[i].xmax - bC[i].x0 - x_0 / non_dimensionalize_coeff), std::abs(bC[i].xmin - bC[i].x0 - x_0 / non_dimensionalize_coeff)), std::min(std::abs(bC[i].ymax - bC[i].y0 - y_0 / non_dimensionalize_coeff), std::abs(bC[i].ymin - bC[i].y0 - y_0 / non_dimensionalize_coeff))) / dx;
	std::string _filename = folder + "/Nu_roundintegtal_x0_" + std::to_string(static_cast<int>(x_0)) + "_y0_" + std::to_string(static_cast<int>(y_0)) + "_" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << y_name << "\t" << "T" << "\t" << "Qizl" << "\t" << "Qel" << "\t" << "Ra_out" << "\t" << "Qck_out" << "\t" << "Alpha" << "\t" << "Nu" << '\n';
#pragma omp parallel for
	for (float j = 1; j < a; j++) {
		float k = 0;
		float Tint = 0;
		float Qizlint = 0;
		float Qelint = 0;
		float Ra_outint = 0;
		float Qck_outint = 0;
		float Alpha_int = 0;
		float Nu_int = 0;
		for (int t = 0; t < T_.size(); t++) {
			if (sqrt((X[t] + x_0 / non_dimensionalize_coeff) * (X[t] + x_0 / non_dimensionalize_coeff) + (Y[t] + y_0 / non_dimensionalize_coeff) * (Y[t] + y_0 / non_dimensionalize_coeff)) <= j*dx / non_dimensionalize_coeff) {
				Tint += T_[t];
				Qizlint += Qizl_[t];
				Qelint += Qel_[t];
				Ra_outint += Ra2_[t];
				Qck_outint += Qck2_[t];
				Alpha_int += Alpha_[t];
				Nu_int += Nu_[t];
				k++;
			}
		}
		if (k > 0) {
			f << j * dx / non_dimensionalize_coeff << "\t" << Tint / k << "\t" << Qizlint / k << "\t" << Qelint / k << "\t" << Ra_outint / k << "\t" << Qck_outint / k << "\t" << Alpha_int / k << "\t" << Nu_int / k << '\n';
		}	
	}
	f.close();
}

