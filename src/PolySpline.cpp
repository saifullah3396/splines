/**
 * @file PolySpline.h
 *
 * This file implements the class PolySpline
 *
 * @author <A href="mailto:saifullah3396@gmail.com">Saifullah</A>
 * @date 12 Aug 2017  
 */

#include "PolySpline.h"

template <typename Scalar>
PolySpline<Scalar>::PolySpline(
  const unsigned& degree,
  const unsigned& dim,
  const Matrix<Scalar, Dynamic, Dynamic>& controlPoints,
  const Matrix<Scalar, Dynamic, 1>& knots, 
  const Scalar& stepSize, 
  const Matrix<Scalar, Dynamic, Dynamic>& boundaryConds) :
  Spline<Scalar>(
    degree,
    dim, 
    controlPoints, 
    knots,
    stepSize,  
    PIECE_WISE_POLY),
  boundaryConds(boundaryConds)
{
}

template <typename Scalar>
PolySpline<Scalar>::PolySpline(
  const unsigned& degree,
  const unsigned& dim,
  const Matrix<Scalar, Dynamic, Dynamic>& controlPoints,
  const Scalar& splineTime, 
  const Scalar& stepSize, 
  const Matrix<Scalar, Dynamic, Dynamic>& boundaryConds) :
  Spline<Scalar>(
    degree,
    dim, 
    controlPoints, 
    splineTime,
    stepSize,  
    PIECE_WISE_POLY),
  boundaryConds(boundaryConds)
{
}

template <typename Scalar>
PolySpline<Scalar>::PolySpline(const string& filePath) : 
  Spline<Scalar>(filePath)
{
  PolySpline<Scalar>::splineFromXml(filePath);
}

template <typename Scalar>
void
PolySpline<Scalar>::validateParameters()
{
  assert(this->controlPoints.rows() == this->knots.size() + 1 && "Control points rows must be of size 'nKnots + 1'");
  assert(this->controlPoints.cols() == this->dim && "Control points cols do not match spline dimensions");
  assert(boundaryConds.rows() == 2 && "Boundary condition rows must be of size '2'");
  assert(boundaryConds.cols() == this->dim && "Boundary condition cols do not match spline dimensions");
}

template <typename Scalar>
void
PolySpline<Scalar>::evaluateSpline(vector<vector<Scalar> >& spline,
  vector<Scalar>& splineTime, const unsigned& derivative)
{
  splineTime.clear();
  Scalar totalTime = 0;
  for (int i = 0; i < this->knots.size(); ++i)
    totalTime += this->knots[i];
  Scalar tTime = 0;
  while (tTime + this->stepSize <= totalTime + 1e-6) {
    tTime += this->stepSize;
    splineTime.push_back(tTime);
  }

  //cout << "totalTime: " << totalTime << endl;
  //cout << "tTime: " << tTime << endl;

  vector<Scalar> times;
  times.resize(this->nKnots + 1);
  times[0] = 0.f;
  for (int i = 1; i < times.size(); ++i) {
    times[i] = times[i - 1] + this->knots[i - 1];
  }

  spline.clear();
  spline.resize(this->dim);
  for (int i = 0; i < splineTime.size(); i++) {
    Scalar t = splineTime[i];
    int knot = 0;
    for (int j = 1; j < times.size(); ++j) {
      if (t <= times[j] + this->stepSize / 2) {
        knot = j - 1;
        break;
      }
    }
    for (int k = 0; k < this->dim; ++k) {
      if (derivative == 0) {
        Scalar pos =
          this->coeffs[0](knot, k) * pow(times[knot + 1] - t, 3) + this->coeffs[1](knot, k) * pow(
            t - times[knot],
            3) + this->coeffs[2](knot, k) * (t - times[knot]) + this->coeffs[3](knot, k) * (times[knot + 1] - t);
        spline[k].push_back(pos);
      } else if (derivative == 1) {
        Scalar vel =
          -3 * this->coeffs[0](knot, k) * pow(times[knot + 1] - t, 2) + 3 * this->coeffs[1](
            knot,
            k) * pow(t - times[knot], 2) + this->coeffs[2](knot, k) - this->coeffs[3](
            knot,
            k);
        spline[k].push_back(vel);
      } else if (derivative == 2) {
        Scalar acc =
          6 * this->coeffs[0](knot, k) * (times[knot + 1] - t) - 6 * this->coeffs[1](
            knot,
            k) * (t - times[knot]);
        spline[k].push_back(acc);
      }
    }
  }

  //cout << "size: " << spline[0].size() << endl;
  //for (int i = 0; i < splineTime.size(); ++i)
  //{
  //	cout << "Time[" << i << "]: " << splineTime[i] << endl;
  //}
}

template <typename Scalar>
void
PolySpline<Scalar>::evaluateCoeffs(const Matrix<Scalar, Dynamic, 1>& knots)
{
  assert(knots.size() == this->nKnots && "Knot sizes do not match");
  this->knots = knots;
  genParams();
}

template <typename Scalar>
void PolySpline<Scalar>::splineFromXml(const string& filePath)
{
  fstream fs;
  fs.open(filePath);
  boost::property_tree::ptree pt;
  read_xml(fs, pt);
  fs.close();
  try {
    BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, pt.get_child("spline")) {
      //! Parsing attributes
      if (v.first == "boundary_conditions") { //! parsing coefficients
        size_t i = 0;
        BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, v.second) {
          if (v.first == "<xmlattr>") {
            boundaryType = v.second.get_child("type").data();
            if (boundaryType != "clamped" && boundaryType != "natural")
              throw SplineException<Scalar>(
                this, 
                "Invalid type specified for boundary conditions.", 
                false);
            boundaryConds.resize(2, this->dim);
          } else {
            if (i >= boundaryConds.rows()) {
              throw SplineException<Scalar>(
                this, 
                "Boundary conditions matrix row size mismatch.", 
                false);
            }
            size_t j = 0;
            BOOST_FOREACH(const boost::property_tree::ptree::value_type &v, v.second) {
              if (j >= this->dim) {
                throw SplineException<Scalar>(
                  this, 
                  "Boundary conditions matrix column size mismatch.", 
                  false
                );
              }
              boundaryConds(i, j) = boost::lexical_cast<Scalar>(v.second.data());
              ++j;
            }
            ++i;
          }
        }
      }
    }
  } catch (SplineException<Scalar> &e) {
    std::cerr << e.what() << std::endl;
  } catch (exception &e) {
    std::cerr << "Failed to parse spline with error: " << e.what() << std::endl;
  }
  validateParameters();
}

template class PolySpline<float>;
template class PolySpline<double>;
