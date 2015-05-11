#pragma once

// Common functions for all PointCloud functions
#include <string.h>
#include <stdio.h>
#include <spii/auto_diff_term.h>
#include <spii/transformations.h>
#include <spii/solver.h>
#include "utils/mexutils.h"

class Projection {
public:

	Projection(unsigned int dimensions) : dimensions(dimensions)
	{}

	// Project point on tangent plane
	template<typename R, typename Q>
	R operator()(const R* const tangent_plane, const Q* const point) const {
		R term = tangent_plane[dimensions];

		for (int i = 0; i < dimensions; ++i) {
			term += point[i]*tangent_plane[i];
		}

		R coeff = tangent_plane[0]*tangent_plane[0];
		for (int i = 1; i < dimensions; ++i) {
			coeff += tangent_plane[i]*tangent_plane[i];
		}

		return term/coeff;
	}

	template<typename R, typename Q>
	R operator()(const R* const tangent_plane, vector<Q> point) const {
		R term = tangent_plane[dimensions];

		for (int i = 0; i < dimensions; ++i) {
			term += point[i]*tangent_plane[i];
		}

		R coeff = tangent_plane[0]*tangent_plane[0];
		for (int i = 1; i < dimensions; ++i) {
			coeff += tangent_plane[i]*tangent_plane[i];
		}

		return term/coeff;
	}

private:
	unsigned int dimensions;
};


class Quadratic_data
{
	public:
		Quadratic_data(unsigned int dimensions, double data_weight, double* point)
		: dimensions(dimensions), data_weight(data_weight), point(point)
		{}

		template<typename R>
 		R operator()(const R* const tangent_plane) const
		{
			R term = tangent_plane[dimensions];
			for (unsigned int i = 0; i < dimensions; ++i) {
				term += point[i]*tangent_plane[i];
			}
			term = term*term;

			R coeff = tangent_plane[0]*tangent_plane[0];
			for (int i = 1; i < dimensions; ++i) {
				coeff += tangent_plane[i]*tangent_plane[i];
			}

			return data_weight*(term/coeff);
		}

protected:
  const unsigned int dimensions;
  const double data_weight;
  const double* const point;
};

class Default_regularization
{
	public:
		Default_regularization(
			unsigned int dimensions, double regularization_weight, double tol, double* point0, double* point1, double eps) :
			dimensions(dimensions), projection(dimensions), regularization_weight(regularization_weight), tol(tol),
			point0(point0), point1(point1), eps(eps)
		{}

		template<typename R>
 		R operator()(const R* const tp0, const R* const tp1) const
		{
			// Project on the tangent planes
			vector<R> p(dimensions);
			vector<R> q(dimensions);

			auto coef1 = projection(tp1, point1);
			auto coef0 = projection(tp0, point0);
			for (unsigned int d = 0; d < dimensions; ++d) {
				p[d] = point0[d] - coef0*tp0[d];
				q[d] = point1[d] - coef1*tp1[d];
			}

			vector<R> dq(dimensions);
			auto coef2 = projection(tp0, q);
			for (unsigned int d = 0; d < dimensions; ++d) {
				dq[d] = coef2*tp0[d];
			}

			R numerator  	=  R(eps);
			R denominator  	=  R(eps);
			for (unsigned int d = 0; d < dimensions; ++d) {
				numerator += (dq[d])*(dq[d]);
				denominator += (p[d] - q[d])*(p[d] - q[d]);
			}

			R cost = numerator / denominator;

			if (cost > R(tol)) {
				cost = tol;
			}

			return regularization_weight*cost;
		}

protected:
  Projection projection;
  const unsigned int dimensions;

  const double regularization_weight;
  const double tol;
  const double eps;

  const double* const point0;
  const double* const point1;
};


class Linear_regularization
{
	public:
		Linear_regularization(
			unsigned int dimensions, double regularization_weight, double tol, double* point0, double* point1, double eps) :
			dimensions(dimensions), projection(dimensions), regularization_weight(regularization_weight), tol(tol), eps(eps),
			point0(point0), point1(point1)
		{}

		template<typename R>
 		R operator()(const R* const tp0, const R* const tp1) const
		{
			// Project on the tangent planes
			vector<R> p(dimensions);
			vector<R> q(dimensions);

			auto coef1 = projection(tp1, point1);
			auto coef0 = projection(tp0, point0);
			for (unsigned int d = 0; d < dimensions; ++d) {
				p[d] = point0[d] - coef0*tp0[d];
				q[d] = point1[d] - coef1*tp1[d];
			}

			vector<R> dq(dimensions);
			auto coef2 = projection(tp0, q);
			for (unsigned int d = 0; d < dimensions; ++d) {
				dq[d] = coef2*tp0[d];
			}

			R numerator  	=  R(eps);
			R denominator  	=  R(eps);
			for (unsigned int d = 0; d < dimensions; ++d) {
				numerator += (dq[d])*(dq[d]);
				denominator += (p[d] - q[d])*(p[d] - q[d]);
			}

			R cost = sqrt(numerator / denominator);

			if (cost > R(tol)) {
				cost = tol;
			}

			return regularization_weight*cost;
		}

protected:
  Projection projection;
  const unsigned int dimensions;

  const double regularization_weight;
  const double tol;
  const double eps;

  const double* const point0;
  const double* const point1;
};


class Qudratic_regularization
{
	public:
		Qudratic_regularization(
			unsigned int dimensions, double regularization_weight, double tol, double* point0, double* point1, double eps) :
			dimensions(dimensions), projection(dimensions), regularization_weight(regularization_weight), tol(tol), eps(eps),
			point0(point0), point1(point1)
		{}

		template<typename R>
 		R operator()(const R* const tp0, const R* const tp1) const
		{
			// Project on the tangent planes
			vector<R> p(dimensions);
			vector<R> q(dimensions);

			auto coef1 = projection(tp1, point1);
			auto coef0 = projection(tp0, point0);
			for (unsigned int d = 0; d < dimensions; ++d) {
				p[d] = point0[d] - coef0*tp0[d];
				q[d] = point1[d] - coef1*tp1[d];
			}

			vector<R> dq(dimensions);
			auto coef2 = projection(tp0, q);
			for (unsigned int d = 0; d < dimensions; ++d) {
				dq[d] = coef2*tp0[d];
			}

			R numerator  	=  R(eps);
			R denominator  	=  R(eps);
			for (unsigned int d = 0; d < dimensions; ++d) {
				numerator += (dq[d])*(dq[d]);
				denominator += (p[d] - q[d])*(p[d] - q[d]);
			}

			R cost = numerator / (denominator*sqrt(denominator));

			if (spii::to_double(cost) > tol) {
				cost = R(tol);
			}

			return regularization_weight*cost;
		}

protected:
  Projection projection;
  const unsigned int dimensions;

  const double regularization_weight;
  const double tol;
  const double eps;

  const double* const point0;
  const double* const point1;
};

class Length_regularization
{
	public:
		Length_regularization(
		 	unsigned int dimensions, double regularization_weight, double tol, double* point0, double* point1, double eps) :
			dimensions(dimensions), projection(dimensions), regularization_weight(regularization_weight), tol(tol), eps(eps),
			point0(point0), point1(point1)
		{}

		template<typename R>
 		R operator()(const R* const tp0, const R* const tp1) const
		{
			// Project on the tangent planes
			vector<R> p(dimensions);
			vector<R> q(dimensions);

			auto coef1 = projection(tp1, point1);
			auto coef0 = projection(tp0, point0);
			for (unsigned int d = 0; d < dimensions; ++d) {
				p[d] = point0[d] - coef0*tp0[d];
				q[d] = point1[d] - coef1*tp1[d];
			}

			R cost  	=  R(eps);
			for (unsigned int d = 0; d < dimensions; ++d) {
				cost += (p[d] - q[d])*(p[d] - q[d]);
			}


			if (spii::to_double(cost) > tol) {
				cost = R(tol);
			}

			return regularization_weight*cost;
		}

protected:
  Projection projection;
  const unsigned int dimensions;

  const double regularization_weight;
  const double tol;
  const double eps;

  const double* const point0;
  const double* const point1;
};