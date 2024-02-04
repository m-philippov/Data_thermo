#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
struct Color {
	float r, g, b;
	Color();
	Color(float r, float g, float b);
	~Color();
};

class Image {
public:
	Image();
	Image(int width, int height);
	~Image();

	Color GetColor(int x, int y) const;
	void SetColor(const Color& color, int x, int y);
	void Read(const std::string& path);
	void Export(const std::string& path) const;
	std::vector<float> ExportTemperature(int& width, int& height) const;
protected:
	int m_width;
	int m_height;
	std::vector<Color> m_colors;
	std::vector<float> temperature;
    
};


/*
	Image image(width, height);
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			image.SetColor(Color(float(x) / float(width), 1.0f - ((float)x / (float)width), (float)y/(float)height),x,y);
		}
	}
	image.Export("image2.bmp");
*/