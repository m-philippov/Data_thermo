#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "../read_begin_conditional.h"
class Effective {
public:
	Effective(const std::vector<float>& T, const float& Tcold, const float& Thot, const std::vector<float>& X, 
		const std::vector<float>& Y, const float& characteristic_size, const std::string& Xname, const std::string& Yname, const std::string& filename);
	void PrintEffective(const std::string& delimetr);


private:
	std::vector<float> _T;
	std::vector<float> _X; // 
	std::vector<float> _Y; //
	std::vector<float> _eff; // эффективность
	float _characteristic_size;
	std::string _X_Name;
	std::string _Y_Name;
	std::string _filename;
};

void AreaEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name);
void AreaIntegralEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name);

void AlongXbEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name);
void AlongXEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float X0);
void AlongYEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float Y0);

void AlongXbIntegralEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float dx);
void AlongXIntegralEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float dx, float y);
void AlongYIntegralEffectiveness(int width_pix, int height_pix, int margin_from_edges_pix, float linear_ax, float linear_b, float x_begin_coord, float y_begin_coord, float non_dimensionalize_coeff, const std::vector<beginConditionals>& bC, std::string x_name, std::string y_name, std::string temp_name, std::string eff_name, float dx, float y);
