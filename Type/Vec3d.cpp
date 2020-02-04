/************************************************************************

	Copyright 2012-2013 Luciano Buglioni

	Contact: luciano.buglioni@gmail.com

	This file is a part of FluxSol

	FluxSol is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    Free CFD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/

#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include "Vec3d.h"
#include "Scalar.h"

using namespace std;

namespace FluxSol
{

//This causes conflict in linking time
//Vec3D::Vec3D(double x=0., double y=0., double z=0.) {
//	comp[0]=x;
//	comp[1]=y;
//	comp[2]=z;
//}

void Vec3D::Set_angles()
{
    double n=fabs(*this);
    for (int i=0;i<3;i++) angle[i]=acos(comp[0]/n);
}

double Vec3D::dot(const Vec3D &right) {
	Vec3D opa;
	opa.comp[0]=comp[0]*right.comp[0];
	return (comp[0]*right.comp[0]+comp[1]*right.comp[1]+comp[2]*right.comp[2]);
}

double fabs(const Vec3D& vec) {
	return sqrt(vec.comp[0]*vec.comp[0]+vec.comp[1]*vec.comp[1]+vec.comp[2]*vec.comp[2]);
}

Vec3D Vec3D::cross(const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]= comp[1]*right.comp[2]-comp[2]*right.comp[1];
	temp.comp[1]=-comp[0]*right.comp[2]+comp[2]*right.comp[0];
	temp.comp[2]= comp[0]*right.comp[1]-comp[1]*right.comp[0];
	return temp;
}


//Equal to inner
const Scalar Vec3D::operator& (const Vec3D &right) const
{
	Scalar temp;
	for (int c=0;c<3;c++) temp.Val()+=comp[c]*right.comp[c];
	return temp;
}

Vec3D operator^ (const Vec3D &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]= left.comp[1]*right.comp[2]-left.comp[2]*right.comp[1];
	temp.comp[1]=-left.comp[0]*right.comp[2]+left.comp[2]*right.comp[0];
	temp.comp[2]= left.comp[0]*right.comp[1]-left.comp[1]*right.comp[0];
	return temp;
}

//Used in Scalar too, and hence in templates
const Vec3D Vec3D::normalize(void) {
	return (*this) /= fabs(*this);
}

const Scalar Vec3D::Norm(void) const {
	//
	//Scalar r = x()*x() + y()*y() + z()*z();
	Vec3D v = *this;
	double r = v.x()*v.x() + v.y()*v.y() + v.z()*v.z();
	return r;
}

//Vec3D& Vec3D::operator= (const Vec3D &right) {
//	comp[0]=right.comp[0];
//	comp[1]=right.comp[1];
//	comp[2]=right.comp[2];
//	return *this;
//}

Vec3D& Vec3D::operator= (const double &right) {
	comp[0]=right;
	comp[1]=right;
	comp[2]=right;
	return *this;
}

Vec3D& Vec3D::operator= (const std::vector<double> &right) {
	comp[0]=right[0];
	comp[1]=right[1];
	comp[3]=right[2];
	return *this;
}

Vec3D &Vec3D::operator*= (const double &right) {
	comp[0]*=right;
	comp[1]*=right;
	comp[2]*=right;
	return *this;
}

const Vec3D Vec3D::operator*(const double &right) const
{
	Vec3D temp;
	temp.comp[0]=comp[0]*right;
	temp.comp[1]=comp[1]*right;
	temp.comp[2]=comp[2]*right;
	return temp;
}

Vec3D &Vec3D::operator/= (const double &right) {
	comp[0]/=right;
	comp[1]/=right;
	comp[2]/=right;
	return *this;
}

const Vec3D Vec3D::operator/ (const double &right) const {
	Vec3D temp;
	temp.comp[0]=comp[0]/right;
	temp.comp[1]=comp[1]/right;
	temp.comp[2]=comp[2]/right;
	return temp;
}

//le sumo un valor constante
Vec3D &Vec3D::operator+= (const double &right) {
	comp[0]+=right;
	comp[1]+=right;
	comp[2]+=right;
	return *this;
}

Vec3D &Vec3D::operator+= (const Vec3D &right) {
	comp[0]+=right.comp[0];
	comp[1]+=right.comp[1];
	comp[2]+=right.comp[2];
	return *this;
}

const Vec3D Vec3D::operator+ (const double &right) const
{
	Vec3D temp;
	temp.comp[0]=comp[0]+right;
	temp.comp[1]=comp[1]+right;
	temp.comp[2]=comp[2]+right;
	return temp;
}


const Vec3D Vec3D::operator+ (const Vec3D &right) const
{
	Vec3D temp;
	temp.comp[0]=comp[0]+right.comp[0];
	temp.comp[1]=comp[1]+right.comp[1];
	temp.comp[2]=comp[2]+right.comp[2];
	return temp;
}

const Vec3D Vec3D::operator- (const Vec3D &right) const
{
	Vec3D temp;
	temp.comp[0]=comp[0]-right.comp[0];
	temp.comp[1]=comp[1]-right.comp[1];
	temp.comp[2]=comp[2]-right.comp[2];
	return temp;
}

Vec3D &Vec3D::operator-= (const double &right) {
	comp[0]-=right;
	comp[1]-=right;
	comp[2]-=right;
	return *this;
}

Vec3D &Vec3D::operator-= (const Vec3D &right) {
	comp[0]-=right.comp[0];
	comp[1]-=right.comp[1];
	comp[2]-=right.comp[2];
	return *this;
}

const Vec3D Vec3D::operator- (const double &right) const
{
	Vec3D temp;
	temp.comp[0]=comp[0]-right;
	temp.comp[1]=comp[1]-right;
	temp.comp[2]=comp[2]-right;
	return temp;
}

Vec3D operator*(const double &left, Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left*right.comp[0];
	temp.comp[1]=left*right.comp[1];
	temp.comp[2]=left*right.comp[2];
	return temp;
}

Vec3D& operator/(const Vec3D &left, Scalar &right) {
	Vec3D temp;
	temp.comp[0]=left.comp[0]/right.Val();
	temp.comp[1]=left.comp[1]/right.Val();
	temp.comp[2]=left.comp[2]/right.Val();
	return temp;
}

const Vec3D operator/ (const double &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left/right.comp[0];
	temp.comp[1]=left/right.comp[1];
	temp.comp[2]=left/right.comp[2];
	return temp;
}

const Vec3D operator+ (const double &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left+right.comp[0];
	temp.comp[1]=left+right.comp[1];
	temp.comp[2]=left+right.comp[2];
	return temp;
}



const Vec3D operator- (const double &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left-right.comp[0];
	temp.comp[1]=left-right.comp[1];
	temp.comp[2]=left-right.comp[2];
	return temp;
}


//typename innerProduct<Vec3D, Vec3D> ::type
//Vec3D::operator& (const Vec3D &right)
//{
//	typename innerProduct<Vec3D, Vec3D> ::type ret;
//
//	//return Vec3D(
//
//}

bool Vec3D::operator== (const Vec3D &right) {
	return (comp[0]==right.comp[0] && comp[1]==right.comp[1] && comp[2]==right.comp[2]);
}

bool Vec3D::operator!= (const Vec3D &right) {
	return (comp[0]!=right.comp[0] || comp[1]!=right.comp[1] || comp[2]!=right.comp[2]);
}


ostream &operator<< (ostream &output,const Vec3D &right) {
	output << "[" << right.comp[0] << "," << right.comp[1] << "," << right.comp[2] << "]";
	return output;
}

//Operadores que implican varias clases



string Vec3D::Imprimir_Coord()
{
    string s,s2;
    char temp[100];
    for (int i=0;i<3;i++){
        sprintf(temp,"%1.20e",comp[i]);s2=temp;
		s+=s2+ "  ";}
    return s;
}
//
//Vec3D  Prom(const vector<Vec3D> &vecs)
//{
//    Vec3D vec(0,0,0);
//    for (int i=0;i<vecs.size();i++)
//        vec+=vecs[i];
//    vec/=vecs.size();
//    return vec;
//}

Vec3D Vec3D::Val()
{
    return (*this);
}

vector<double> Vec3D::Comp()
{
	vector<double> r;
	for (int i=0;i<3;i++)
	{
		r.push_back(this->comp[i]);
	}
	return r;
}

const Vec3D operator* (const Scalar &sc, const Vec3D &vector)
{
	Vec3D v;
	Scalar s=sc;
	for (int c=0;c<3;c++)
		v.comp[c]=s.Val()*vector.comp[c];

	return v;
}

const Vec3D operator* (const Vec3D &vector, const Scalar &sc)
{
	Vec3D v;
	Scalar s = sc;
	for (int c=0;c<3;c++)
		v.comp[c]=s.Val()*vector.comp[c];

	return v;
}

Vec3D operator/ (Vec3D &vector, Scalar &sc)
{
	Vec3D v;
	for (int c=0;c<3;c++)
		v.comp[c]=vector.comp[c]/sc.Val();

	return v;
}


Scalar operator& (Vec3D &left, Vec3D &right)
{

	double d=0.;
	for (int c=0;c<3;c++)
		d+=left.comp[c]*right.comp[c];
	Scalar s(d);
	return s;
}

const Vec3D operator* (const double &left, const Vec3D &right)
{
	Vec3D r;
	for (int c=0;c<3;c++)
		r.comp[c]=left*right.comp[c];
	return r;
}

}
