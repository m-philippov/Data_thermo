#pragma once
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
//number expiriment	t1	t2	t3	t4	u	g	m	x0	y0	x1	y1	xmin	xmax	ymin	ymax dx

struct beginConditionals {
	std::string number;
	float t1;
	float t2;
	float t3;
	float t4;
	float u;
	float g;
	float m;
	int x0;
	int y0;
	int x1;
	int y1;
	int xb;
	int xmax;
	int ymin;
	int ymax;
	float dx_mm;
};

std::vector<beginConditionals> InputBeginConditionals(const std::string& nameBeginConditionalFile);
/*std::ostream& operator<<(std::ostream& os, const beginConditionals& bC)
{
	return os << "number= " << bC.number << " t1 = " << bC.t1 << " t2 = " << bC.t2 << " t3 = " << bC.t3 << " t4 = " << bC.t4
		<< " u = " << bC.u << " g = " << bC.g << " m = " << bC.m << " x0 = " << bC.x0 << " y0 = " << bC.y0
		<< " x1 = " << bC.x1 << " y1 = " << bC.y1 << " xmin = " << bC.xmin << " xmax = " << bC.xmax << " ymin = " << bC.ymin
		<< " ymax = " << bC.ymax << " dx = " << bC.dx_mm;
}*/
//for impiging jets Nu
//number expiriment	numjets	s/d	h/d	Q,g/s	U	A	W	tair	tjet	x0	y0	ymin	ymax dx 

struct beginConditionalsNu4 {
	std::string number;
	int num_jets;
	float s_d;
	float h_d;
	float Qg_c;
	float U;
	float A;
	float W;
	float Tair;
	float Tjet;
	int x0;
	int y0;
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	float dx_mm;
};

std::vector<beginConditionalsNu4> InputBeginConditionalsNu4(const std::string& nameBeginConditionalFile);

struct beginConditionalsNu3 {
	std::string number;
	int num_jets;
	float s_d;
	float h_d;
	float Qg_c;
	float U;
	float A;
	float W;
	float Tair;
	float Tjet;
	int x0;
	int y0;
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	float dx_mm;
	float velosity;
};

std::vector<beginConditionalsNu3> InputBeginConditionalsNu3(const std::string& nameBeginConditionalFile);