#pragma once
#include <map>
#include <vector>
enum class Figure : int {
	LINE = 0,
	ROUND,
	RECTANGULAR,
	ARC, //дуга
};
struct MapFigure{
	Figure type;
	std::pair<float, float> point1;
	std::pair<float, float> point2;
	std::pair<float, float> point3;
};
//std::pair<float, float> point_center;
//std::vector<MapFigure> figure_difficult;
