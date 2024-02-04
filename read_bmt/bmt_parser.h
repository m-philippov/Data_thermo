#pragma once
#include "bmt_extractor.h"


class DataExtract : public Image {
public:
	DataExtract(const std::string& a);
	std::vector<Color> DataExportImage();
	std::vector<float> DataExportTemperature();
	void DataPrintGraphTable() const;
	void DataPrintImage() const;
private:
	std::string file_name;
};
