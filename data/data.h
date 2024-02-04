#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cmath>

#define PROPERTIES_OF_AIR(func) OpenPropertiesOFAir(func, #func)
//��������� ���������
const float g_free = 9.8;
const float T0 = 273.15;
const float sigma = 5.67 * static_cast<float>(pow(10.0, -8));


//��������� ���� � ������ ������ � ������ ���
std::map<float, float> AirCriterium(std::ifstream&);

//������������
std::map<float, float> Linealize(const std::map<float, float>&);


void OpenPropertiesOFAir(std::map<float, float>&, const std::string&);

//���������� ������ "���������������� ������� � ����������� �� ����������� � �������"
//AirLambdaPatametr


//���������� ������ "��������� ������� � ����������� �� ����������� � �������"
// AirDensity

//���������� ������ "������������  �������� ������� � ����������� �� ����������� � �������"
// AirDynamicViscosity

//���������� ������ "��������������  �������� ������� � ����������� �� ����������� � �������"
//AirKinematicViscosity

//���������� ������ "�������� ������������ ������� ��� ��������� ������������ � �������"
//AirSpecificHeat

//���������� ������ "��������� ������� � ����������� �� ����������� � �������"
//AirThermalConductivity
// 
//���������� ������ "������������  �������� ������� � ����������� �� ����������� � �������"
//AirThermalDiffusivity

//���������� ������ "��������������  �������� ������� � ����������� �� ����������� � �������"
//AirVolumetric�oefficientofThermalExpansion

//���������� ������ "�������� ������������ ������� ��� ��������� ������������ � �������"
//AirPrandtl

/**/
