#pragma once
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
//номер эксперимента	t1	t2	t3	t4	u	g	m	x0	y0	x1	y1	xmin	xmax	ymin	ymax dx

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
