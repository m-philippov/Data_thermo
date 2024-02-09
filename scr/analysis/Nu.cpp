#include "Nu.h"


//

//	PROPERTIES_OF_AIR(AirLambdaPatametr);
//	PROPERTIES_OF_AIR(AirKinematicViscosity);
//	PROPERTIES_OF_AIR(AirDensity);
//	PROPERTIES_OF_AIR(AirSpecificHeat);
//	PROPERTIES_OF_AIR(AirThermalDiffusivity);
//	PROPERTIES_OF_AIR(AirPrandtl);
//	PROPERTIES_OF_AIR(AirLambdaPatametr);
//	for (int i = 0; i < T.size(); ++i) {
//		_T[i] = T[i];
//		_Tin[i] = (Tair + N * (X_begin - X[i])) / (LinearInterpolation(AirSpecificHeat, T[i]) *
//			LinearInterpolation(AirDensity, T[i]) * velocity * high * lenght * weight);
//		if (_Tin[i] >= T[i]) {
//			_Tin[i] = T[i];
//		}
//		_Qel[i] = N / (_NumberofWall * weight * lenght);
//		_Qizl[i] = sigma * epsilon_out * (std::pow((_T[i] + T0), 4) - std::pow((Tair + T0), 4)) +
//			sigma * epsilon_in * (std::pow((_T[i] + T0), 4) - std::pow((Tair + T0), 4));
//		float Tin0 = _T[i] - 1 / 3 * _Tin[i];
//		float Tin1 = _T[i] - 1 / 3 * Tair;
//		_Ra_in[i] = 2.0 * g_free * (T[i] - _Tin[i]) * pow((X[i]), 3) / (_T[i] + _Tin[i] + 2.0 * T0) / 
//			LinearInterpolation(AirKinematicViscosity, Tin0) / LinearInterpolation(AirThermalDiffusivity, Tin0);
//		_Ra_out[i] = 2.0 * g_free * (T[i] - Tair) * pow((X[i]), 3) / (_T[i] + Tair + 2.0 * T0) / 
//			LinearInterpolation(AirKinematicViscosity, Tin1) / LinearInterpolation(AirThermalDiffusivity, Tin1);
//		_Qck_in[i] = 0.508 * LinearInterpolation(AirLambdaPatametr, Tin0) * (T[i] - _Tin[i]) * 
//			pow((_Ra_in[i] * LinearInterpolation(AirPrandtl, Tin0) / (0.952 + LinearInterpolation(AirPrandtl, Tin0))), 0.25) / X[i];
//		_Qck_out[i] = 0.508 * LinearInterpolation(AirLambdaPatametr, Tin1) * (T[i] - Tair) *
//			pow((_Ra_in[i] * LinearInterpolation(AirPrandtl, Tin1) / (0.952 + LinearInterpolation(AirPrandtl, Tin1))), 0.25) / X[i];
//		_Qw[i] = _Qel[i] - _Qizl[i] - _Qck_in[i] - _Qck_out[i];
//		_Alpha[i] = _Qw[i] / (T[i] - _Tin[i] + 0.001);
//		_Nu[i] = diameter * _Alpha[i] / LinearInterpolation(AirLambdaPatametr, Tin0);
//		if (_Tin[i] >= T[i]) {
//			_Alpha[i] = 0;
//			_Nu[i] = 0;
//		}
//	}
//	
//}
//

Nu::Nu(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionalsNu4>& bC, const std::string& x_name, const std::string& y_name, int i, float lenght_plate, float weight_plate, float epsilon_kraska, float epsilon_material, float diameter_jet)
: width_pix(width_pix), height_pix(height_pix), margin_from_edges_pix(margin_from_edges_pix), linear_ax(linear_ax), linear_b(linear_b), y_begin_coord(y_begin_coord), x_begin_coord(x_begin_coord),non_dimensionalize_coeff(non_dimensionalize_coeff), bC(bC), x_name(x_name), y_name(y_name), i(i)
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
    std::ifstream infile(bC[i].number + ".txt");
    for (int j = 0; j < height_pix * width_pix; j++) {
        infile >> Ti;
		if ((j % width_pix + margin_from_edges_pix < bC[i].xmax) && (j % width_pix - margin_from_edges_pix > bC[i].xmin) && (j / width_pix - margin_from_edges_pix < bC[i].ymax) && (j / width_pix + margin_from_edges_pix > bC[i].ymin)) {
#pragma omp parallel
			{
				T_.push_back(linear_ax * Ti - linear_b);
				float AirLambdaTk = LinearInterpolation(AirLambdaPatametr, T_[k]);
				float AirPrandtlTk = LinearInterpolation(AirPrandtl, T_[k]);
				X.push_back(((j % width_pix) * bC[i].dx_mm - bC[i].x0 * bC[i].dx_mm + x_begin_coord) / non_dimensionalize_coeff);
				Y.push_back((bC[i].y0 * bC[i].dx_mm - (j / width_pix) * bC[i].dx_mm + y_begin_coord) / non_dimensionalize_coeff);
				Qizl_.push_back(sigma * (epsilon_kraska + epsilon_material) * (pow((T_[k] + 273.15), 4) - pow((273.15 + bC[i].Tair), 4)));
				Qel_.push_back(bC[i].W / (weight_plate * lenght_plate));
				//Ra_.push_back(LinearInterpolation(AirVolumetric—oefficientofThermalExpansion, T_[k]) * g_free * (Qel_[k] - Qizl_[k]) * pow(diameter_jet, 4) / AirLambdaTk / LinearInterpolation(AirKinematicViscosity, T_[k]) / LinearInterpolation(AirThermalDiffusivity, T_[k]));
				Ra2_.push_back(AirVolumetric—oefficientofThermalExpansionTair * g_free * (Qel_[k] - Qizl_[k]) * pow(diameter_jet, 4) / AirLambdaTair / AirKinematicViscosityTair / AirThermalDiffusivityTair);
				//Qck1_.push_back(0.508 * AirLambdaTk * (T_[k] - bC[i].Tair) * pow((Ra_[k] * AirPrandtlTk / (0.952 + AirPrandtlTk)), 0.25) / diameter_jet);
				Qck2_.push_back(0.508 * AirLambdaTair * (T_[k] - bC[i].Tair) * pow((Ra2_[k] * AirPrandtlTair / (0.952 + AirPrandtlTair)), 0.25) / diameter_jet);
				Qw_.push_back((Qel_[k] - Qizl_[k] - Qck2_[k]) >= 0 ? (Qel_[k] - Qizl_[k] - Qck2_[k]) : 0);
				Nu_.push_back(diameter_jet * (Qw_[k]) / (T_[k] - bC[i].Tjet) / AirLambdaTk);
				Alpha_.push_back(Qw_[k] / (T_[k] - bC[i].Tjet));
				k++;
			}
        }
    }
    infile.close();
}

void Nu::AreaNu()
{
	std::string _filename = "Nu_area_" + bC[i].number + ".txt";

	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << x_name << ", " << y_name << ", " << "T" << ", " << "Qizl" << ", " << "Qel" << ", " << "Ra_out" << ", " << "Qck_out" << ", " << "Alpha" << ", " << "Nu" << '\n';
#pragma omp parallel for
	for (int j = 0; j < T_.size(); j++) {
		f << X[j] << ", " << Y[j] << ", " << T_[j] << ", " << Qizl_[j] << ", " << Qel_[j] <<  ", " << Ra2_[j] << ", "  << Qck2_[j] << ", " << Alpha_[j] << ", " << Nu_[j] << '\n';
	}
	f.close();
}

void Nu::AreaIntegralNu()
{
	std::string _filename = "Nu_area_integral" + bC[i].number + ".txt";
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
	std::string _filename = "Nu_alongXb_" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << x_name << ", " << y_name << ", " << "T" << ", " << "Qizl" << ", " << "Qel" << ", " << "Ra_in" << ", " << "Ra_out" << ", " << "Qck_in" << ", " << "Qck_out" << ", " << "Alpha" << ", " << "Nu" << '\n';
#pragma omp parallel for
	for (int j = 0; j < T_.size(); j++) {
		if ((X[j] >= bC[i].dx_mm / non_dimensionalize_coeff/2) && (X[j] <= -bC[i].dx_mm / non_dimensionalize_coeff / 2)) {
			f << X[j] << ", " << Y[j] << ", " << T_[j] << ", " << Qizl_[j] << ", " << Qel_[j] << ", " << Ra_[j] << ", " << Ra2_[j] << ", " << Qck1_[j] << ", " << Qck2_[j] << ", " << Alpha_[i] << ", " << Nu_[i] << '\n';
		}
	}
	f.close();
}

void Nu::AlongXNu(float x_0)
{
	std::string _filename = "Nu_alongX_" + std::to_string(static_cast<int>(x_0)) + "mm" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << x_name << ", " << y_name << ", " << "T" << ", " << "Qizl" << ", " << "Qel" << ", " << "Ra_in" << ", " << "Ra_out" << ", " << "Qck_in" << ", " << "Qck_out" << ", " << "Alpha" << ", " << "Nu" << '\n';
#pragma omp parallel for
	for (int j = 0; j < T_.size(); j++) {
		if ((X[j] - x_0 >= bC[i].dx_mm / non_dimensionalize_coeff / 2) && (X[j] - x_0 <= -bC[i].dx_mm / non_dimensionalize_coeff / 2)) {
			f << X[j] << ", " << Y[j] << ", " << T_[j] << ", " << Qizl_[j] << ", " << Qel_[j] << ", " << Ra_[j] << ", " << Ra2_[j] << ", " << Qck1_[j] << ", " << Qck2_[j] << ", " << Alpha_[i] << ", " << Nu_[i] << '\n';
		}
	}
	f.close();
}

void Nu::AlongYNu(float y_0)
{
	std::string _filename = "Nu_alongX_" + std::to_string(static_cast<int>(y_0)) + "mm" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << x_name << ", " << y_name << ", " << "T" << ", " << "Qizl" << ", " << "Qel" << ", " << "Ra_in" << ", " << "Ra_out" << ", " << "Qck_in" << ", " << "Qck_out" << ", " << "Alpha" << ", " << "Nu" << '\n';
#pragma omp parallel for
	for (int j = 0; j < T_.size(); j++) {
		if ((Y[j] - y_0 >= bC[i].dx_mm / non_dimensionalize_coeff / 2) && (Y[j] - y_0 <= -bC[i].dx_mm / non_dimensionalize_coeff / 2)) {
			f << X[j] << ", " << Y[j] << ", " << T_[j] << ", " << Qizl_[j] << ", " << Qel_[j] << ", " << Ra_[j] << ", " << Ra2_[j] << ", " << Qck1_[j] << ", " << Qck2_[j] << ", " << Alpha_[i] << ", " << Nu_[i] << '\n';
		}
	}
	f.close();
}

void Nu::AlongXbIntegralNu(float dx)
{
	float Xint = 0;
	float Yint = 0;
	float Tint;
	float Qizlint;
	float Qelint;
	float Ra_inint;
	float Ra_outint;
	float Qck_inint;
	float Qck_outint;
	float Alpha_int;
	float Nu_int;
	std::string _filename = "Nu_alongXb_" + bC[i].number + ".txt";
	std::ofstream f;
	f.open(_filename, std::ios::out);
	f << x_name << ", " << y_name << ", " << "T" << ", " << "Qizl" << ", " << "Qel" << ", " << "Ra_in" << ", " << "Ra_out" << ", " << "Qck_in" << ", " << "Qck_out" << ", " << "Alpha" << ", " << "Nu" << '\n';
#pragma omp parallel for
	int k = 0;
	for (int j = 0; j < T_.size(); j++) {
		if ((X[j] >= bC[i].dx_mm / non_dimensionalize_coeff / 2 - dx) && (X[j] <= -bC[i].dx_mm / non_dimensionalize_coeff / 2 + dx)) {
			Xint += X[j];
			Yint += Y[j];
			f << X[j] << ", " << Y[j] << ", " << T_[j] << ", " << Qizl_[j] << ", " << Qel_[j] << ", " << Ra_[j] << ", " << Ra2_[j] << ", " << Qck1_[j] << ", " << Qck2_[j] << ", " << Alpha_[i] << ", " << Nu_[i] << '\n';
		}
	}
	f.close();
}

void Nu::AlongXIntegralNu(float dx, float x_0)
{
}

void Nu::AlongYIntegralNu(float dx, float y_0)
{
}

void Nu::RoundIntegralNu(float x_0, float y_0, float dx)
{
}

