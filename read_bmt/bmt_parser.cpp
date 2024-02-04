#include "bmt_parser.h"

DataExtract::DataExtract(const std::string& a) : file_name(a)
{
	Read(file_name);
}

std::vector<Color> DataExtract::DataExportImage()
{
	return m_colors;
}

std::vector<float> DataExtract::DataExportTemperature()
{
	return ExportTemperature(m_width, m_height);
}

void DataExtract::DataPrintGraphTable() const
{
	
	std::string path;
	for (int i = 0; i < file_name.size() - 4; ++i) {
		path += file_name[i];
	}
	path += "tabletemperature.txt";
	std::ofstream f;
	f.open(path, std::ios::out);

	for (int i = 0; i < m_height; i++) {
		for (int j = 0; j < m_width; j++) {
			f << temperature[i*m_width + j] << '\t';
		}
		f << '\n';
	}
	f.close();
}

void DataExtract::DataPrintImage() const
{
	std::string path;
	for (int i = 0; i < file_name.size() - 4; ++i) {
		path += file_name[i];
	}
	path += "image.bmp";
	Export(path);
}