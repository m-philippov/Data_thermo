//#include "Nu.h"
//
////Nu_1Builder::operator Nu_1() const
////{
//	// throw exception if not valid
////	return { _Tair, _epsilon_in, _epsilon_out, _diameter,
////	_Tbegin, _N, _T, _X, _Y, _characteristic_size, _X_Name, _Y_Name, _X_begin,
////	_high, _weight, _lenght, _velocity, _filename, _NumberofWall };
////}
//
//Nu_1Builder& Nu_1Builder::SetTair(float Tair) {
//	_Tair = Tair;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetEpsilonIn(float epsilon_in) {
//	_epsilon_in = epsilon_in;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetEpsilonOut(float epsilon_out) {
//	_epsilon_out = epsilon_out;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetDiameter(float diameter) {
//	_diameter = diameter;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetTbegin(float Tbegin) {
//	_Tbegin = Tbegin;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetN(float N) {
//	_N = N;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetCharacteristicSize(float characteristic_size) {
//	_characteristic_size = characteristic_size;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetXName(float X_Name) {
//	_X_Name = X_Name;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetYName(float Y_Name) {
//	_Y_Name = Y_Name;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetXbegin(float X_begin) {
//	_X_begin = X_begin;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetHigh(float high) {
//	_high = high;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetWeight(float weight) {
//	_weight = weight;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetLenght(float lenght) {
//	_lenght = lenght;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetVelocity(float velocity) {
//	_velocity = velocity;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetFilename(float filename) {
//	_filename = filename;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetNumberOFWall(float NumberofWall) {
//	_NumberofWall = NumberofWall;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetT(std::vector<float>& T) {
//	_T = T;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetX(std::vector<float>& X) {
//	_X = X;
//	return *this;
//}
//Nu_1Builder& Nu_1Builder::SetY(std::vector<float>& Y) {
//	_Y = Y;
//	return *this;
//}
//Nu_1::Nu_1(float Tair, float epsilon_in, float epsilon_out, float diameter, float Tbegin, float N,
//	const std::vector<float>& T, const std::vector<float>& X, const std::vector<float>& Y, float characteristic_size,
//	std::string X_Name, std::string Y_Name, float X_begin, float high, float weight, float lenght,
//	float velocity, std::string filename, int NumberofWall)
//{
//	for (int i = 0; i < filename.size() - 4; ++i) {
//		_filename += filename[i];
//	}
//	_characteristic_size = characteristic_size;
//	_filename += "Nu.txt";
//	_X.resize(X.size());
//	for (int i = 0; i < X.size(); i++) {
//		_X[i] = X[i] / characteristic_size;
//	}
//	_Y.resize(Y.size());
//	for (int i = 0; i < Y.size(); i++) {
//		_T[i] = Y[i] / characteristic_size;
//	}
//	_NumberofWall = NumberofWall;
//	_X_Name = X_Name;
//	_Y_Name = Y_Name;
//	_T.resize(T.size());
//	for (int i = 0; i < T.size(); i++) {
//		_T[i] = T[i];
//	}
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
//void Nu_1::PrintNu_1(const std::string& delimetr)
//{
//	std::ofstream f;
//	f.open(_filename, std::ios::out);
//	f << _X_Name << delimetr << _Y_Name << delimetr << "T" << delimetr << "Tin" << delimetr << "Qel" << delimetr << "Qizl" << delimetr << 
//		"Ra_in" << delimetr << "Qck_in" << delimetr << "Ra_out" << delimetr << "Qck_out" << delimetr << "Alpha" << delimetr << "Nu" << '\n';
//	for (int i = 0; i < _T.size(); i++) {
//		f << _X[i] << delimetr << _Y[i] << delimetr << _T[i] << delimetr << _Tin[i] << delimetr << _Qel[i] << delimetr << _Qizl[i] << delimetr <<
//			_Ra_in[i] << delimetr << _Qck_in[i] << delimetr << _Ra_out[i] << delimetr << _Qck_out[i] << delimetr << _Alpha[i] << delimetr << _Nu[i] << '\n';
//	}
//	f.close();
//}
//
