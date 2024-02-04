#include "bmt_extractor.h"

Color::Color() :r(0), g(0), b(0)
{
}

Color::Color(float r, float g, float b) :r(r), g(g), b(b)
{

}

Color::~Color()
{
}

Image::Image() : m_width(0), m_height(0), m_colors(std::vector<Color>(0))
{
}

Image::Image(int width, int height) : m_width(width), m_height(height), m_colors(std::vector<Color>(width* height))
{

}

Image::~Image()
{
}

Color Image::GetColor(int x, int y) const
{
	return m_colors[y * m_width + x];
}

void Image::SetColor(const Color& color, int x, int y)
{
	m_colors[y * m_width + x].r = color.r;
	m_colors[y * m_width + x].g = color.g;
	m_colors[y * m_width + x].b = color.b;
}

void Image::Read(const std::string& path)
{
	temperature.clear();
	std::ifstream f;
	f.open(path, std::ios::in | std::ios::binary);
	if (!f.is_open()) {
		std::cout << "File could not be opened\n";
		return;
	}
	const int fileHeaderSize = 14;
	const int informationHeaderSize = 40;
	unsigned char fileHeader[fileHeaderSize];
	f.read(reinterpret_cast<char*>(fileHeader), fileHeaderSize);
	if (fileHeader[0] != 'B' || fileHeader[1] != 'M') {
		std::cout << "This specified is not a bitmap image" << std::endl;
		f.close();
		return;
	}
	unsigned char informationHeader[informationHeaderSize];
	f.read(reinterpret_cast<char*>(informationHeader), informationHeaderSize);

	int fileSize = fileHeader[2] + (fileHeader[3] << 8) + (fileHeader[4] << 16) + (fileHeader[5] << 24);
	std::cout << "sileSize = " << fileSize << std::endl;
	m_width = informationHeader[4] + (informationHeader[5] << 8)+(informationHeader[6] << 16) + (informationHeader[7] << 24);
	m_height = informationHeader[8] + (informationHeader[9] << 8) +(informationHeader[10] << 16) + (informationHeader[11] << 24);
	std::cout << "width = " << m_width << std::endl;
	std::cout << "height = " << m_height << std::endl;
	this->temperature.resize(m_width* m_height);
	m_colors.resize(m_width * m_height);
	const int paddingAmount = ((4 - (m_width * 3) % 4) % 4);
	int number_char = 54;
	for (int y = 0; y < m_height; y++) {
		for (int x = 0; x < m_width; x++) {
			unsigned char color[3];
			f.read(reinterpret_cast<char*>(color), 3);
			number_char += 3;
			m_colors[y * m_width + x].r = static_cast<float>(color[2]) / 255.0f;
			m_colors[y * m_width + x].g = static_cast<float>(color[1]) / 255.0f;
			m_colors[y * m_width + x].b = static_cast<float>(color[0]) / 255.0f;
		}
		f.ignore(paddingAmount);
	}
	unsigned char fileReader_1[1];
	int leveling = number_char % 16;
	for (int i = 0; i < 16 - leveling - 1; i++) {
		f.read(reinterpret_cast<char*>(fileReader_1), 1);
		number_char++;
	}
	//leveling = number_char % 16;
	unsigned char fileReader_4[4];
	int a, b;
	//if ()
	int c = 1000000 + m_width * m_height;
	unsigned char* fileReader_c = new unsigned char[c];
	f.read(reinterpret_cast<char*>(fileReader_c), c);
	number_char += c;
	delete[] fileReader_c;
	do {
		do {
			f.read(reinterpret_cast<char*>(fileReader_4), 4);
			a = fileReader_4[0] + (fileReader_4[1] << 8) + (fileReader_4[2] << 16) + (fileReader_4[3] << 24);
			
			number_char += 4;
		} while ((m_width != a) );
		f.read(reinterpret_cast<char*>(fileReader_4), 4);
		b = fileReader_4[0] + (fileReader_4[1] << 8) + (fileReader_4[2] << 16) + (fileReader_4[3] << 24);
		number_char += 4;
	} while ((m_height != b));
	
	
	unsigned int p1;
	float p;
	
	for (int i = 0; i < m_height * m_width; i++) {
		f.read(reinterpret_cast<char*>(fileReader_4), 4);
		p1 = fileReader_4[0] + (fileReader_4[1] << 8) + (fileReader_4[2] << 16) + (fileReader_4[3] << 24);
		
		temperature[i] =  *((float*)&p1);
	}
	f.close();
	std::cout << "File read\n";
}

void Image::Export(const std::string& path) const
{
	std::ofstream f;
	f.open(path, std::ios::out | std::ios::binary);

	if (!f.is_open()) {
		std::cout << "File could not be opened\n";
		return;
	}

	unsigned char bmpPad[3] = { 0,0,0 };
	const int paddingAmount = ((4 - (m_width * 3) % 4) % 4);
	const int fileHeaderSize = 14;
	const int informationHeaderSize = 40;
	const int fileSize = fileHeaderSize + informationHeaderSize + m_width * m_height * 3 + paddingAmount * m_width;
	unsigned char fileHeader[fileHeaderSize];
	// File type
	fileHeader[0] = 'B';
	fileHeader[1] = 'M';
	// File size
	fileHeader[2] = fileSize;
	fileHeader[3] = fileSize >> 8;
	fileHeader[4] = fileSize >> 16;
	fileHeader[5] = fileSize >> 24;
	// Reserved 1 (Not used)
	fileHeader[6] = 0;
	fileHeader[7] = 0;
	// Reserved 2 (Not used)
	fileHeader[8] = 0;
	fileHeader[9] = 0;
	// Pixel data offset
	fileHeader[10] = fileHeaderSize + informationHeaderSize;
	fileHeader[11] = 0;
	fileHeader[12] = 0;
	fileHeader[13] = 0;
	unsigned char informationHeader[informationHeaderSize];
	// Header size
	informationHeader[0] = informationHeaderSize;
	informationHeader[1] = 0;
	informationHeader[2] = 0;
	informationHeader[3] = 0;
	// Image width
	informationHeader[4] = m_width;
	informationHeader[5] = m_width >> 8;
	informationHeader[6] = m_width >> 16;
	informationHeader[7] = m_width >> 24;
	// Image height
	informationHeader[8] = m_height;
	informationHeader[9] = m_height >> 8;
	informationHeader[10] = m_height >> 16;
	informationHeader[11] = m_height >> 24;
	// Planes
	informationHeader[12] = 1;
	informationHeader[13] = 0;
	// Bits per pixel (RGB)
	informationHeader[14] = 24;
	informationHeader[15] = 0;
	// Compression (No compression)
	informationHeader[16] = 0;
	informationHeader[17] = 0;
	informationHeader[18] = 0;
	informationHeader[19] = 0;
	// Image size (No compression)
	informationHeader[20] = 0;
	informationHeader[21] =	0;
	informationHeader[22] = 0;
	informationHeader[23] = 0;
	// X pixels per meter (Not specified)
	informationHeader[24] = 0;
	informationHeader[25] = 0;
	informationHeader[26] =
	informationHeader[27] = 0;
	// Y pixels per meter (Not specified)
	informationHeader[28] = 0;
	informationHeader[29] = 0;
	// X pixels per meter (Not specified)
	informationHeader[24] = 0;
	informationHeader[25] = 0;
	informationHeader[26] = 0;
	informationHeader[27] = 0;
	//	Y pixels per meter(Not specified)
	informationHeader[28] = 0;
	informationHeader[29] = 0;
	informationHeader[30] = 0;
	informationHeader[31] = 0;
	// Total colors (Color palette not used)
	informationHeader[32] =	0;
	informationHeader[33] = 0;
	informationHeader[34] = 0;
	informationHeader[35] = 0;
	// Important colors (Generally ignored)
	informationHeader[36] =	0;
	informationHeader[37] = 0;
	informationHeader[38] = 0;
	informationHeader[39] = 0;
	f.write(reinterpret_cast<char*>(fileHeader), fileHeaderSize);
	f.write(reinterpret_cast<char*>(informationHeader), informationHeaderSize);
	for (int y = 0; y < m_height; y++) {
		for (int x = 0; x < m_width; x++) {
			unsigned char r = static_cast<unsigned char>(GetColor(x, y).r * 255.0f);
			unsigned char g = static_cast<unsigned char>(GetColor(x, y).g * 255.0f);
			unsigned char b = static_cast<unsigned char>(GetColor(x, y).b * 255.0f);
			unsigned char color[] = { b, g, r };
			f.write(reinterpret_cast<char*>(color), 3);
		}
		f.write(reinterpret_cast<char*>(bmpPad), paddingAmount);
	}
	f.close();
	std::cout << "File created\n";
}

std::vector<float> Image::ExportTemperature(int& width, int& height) const
{
	width = m_width;
	height = m_height;
	return temperature;
}
