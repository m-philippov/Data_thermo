#include "interpolate.h"

float LinearInterpolation(std::map<float, float> a, float temperature)
{
	{
		float k = 10 * std::floor(temperature / 10);
		auto f = a.find(k);
		if (f != a.end())
			return a[k] + (a[10 * std::ceil(temperature / 10)] - a[k]) * (temperature / 10 - k / 10);
		else {
			throw std::invalid_argument("outside the temperature table range");
		}
	}
}
