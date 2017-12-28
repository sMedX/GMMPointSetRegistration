#ifndef itkGMML2PointSetToPointSetMetric_hxx
#define itkGMML2PointSetToPointSetMetric_hxx

#include "itkGMML2PointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMML2PointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMML2PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();

  this->m_NormalizeFactor = 1.0 / this->m_MovingPointSet->GetNumberOfPoints();
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointType & point) const
{
  const double factor1 = this->m_TransformedPointSet->GetNumberOfPoints() * this->m_FixedPointSet->GetNumberOfPoints();
  const double factor2 = this->m_TransformedPointSet->GetNumberOfPoints() * this->m_TransformedPointSet->GetNumberOfPoints();

  const double scale1 = 0.5 * (this->m_FixedPointSetScale * this->m_FixedPointSetScale + this->m_MovingPointSetScale * this->m_MovingPointSetScale);
  const double scale2 = this->m_MovingPointSetScale * this->m_MovingPointSetScale;

  // compute value for the first sum
  double value1 = 0;
 
  for (FixedPointIterator it = this->m_FixedPointSet->GetPoints()->Begin(); it != this->m_FixedPointSet->GetPoints()->End(); ++it) {
    const double distance = point.SquaredEuclideanDistanceTo(it.Value());
    const double expval = std::exp(-distance / scale1);
    value1 += expval;
  }

  // compute value for the second sum
  double value2 = 0;
 
  for (MovingPointIterator it = this->m_TransformedPointSet->GetPoints()->Begin(); it != this->m_TransformedPointSet->GetPoints()->End(); ++it) {
    const double distance = point.SquaredEuclideanDistanceTo(it.Value());
    const double expval = std::exp(-distance / scale2);
    value2 += expval;
  }

  // local value
  const double value = value2 / factor2 - 2.0 * value1 / factor1;

  return value;
}

template<typename TFixedPointSet, typename TMovingPointSet>
void
GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative(const MovingPointType & point, MeasureType & value, LocalDerivativeType & derivative) const
{
  const double factor1 = this->m_FixedPointSet->GetNumberOfPoints();
  const double factor2 = this->m_TransformedPointSet->GetNumberOfPoints();

  const double scale1 = 0.5 * (this->m_FixedPointSetScale * this->m_FixedPointSetScale + this->m_MovingPointSetScale * this->m_MovingPointSetScale);
  const double scale2 = this->m_MovingPointSetScale * this->m_MovingPointSetScale;

  // compute value and derivative gradient for the first sum
  double value1 = 0;
  LocalDerivativeType derivative1;
  derivative1.Fill(NumericTraits<LocalDerivativeValueType>::ZeroValue());

  for (FixedPointIterator it = this->m_FixedPointSet->GetPoints()->Begin(); it != this->m_FixedPointSet->GetPoints()->End(); ++it) {
    const typename FixedPointSetType::PointType fixedPoint = it.Value();
    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / scale1);
    value1 += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) {
      derivative1[dim] += expval * (point[dim] - fixedPoint[dim]);
    }
  }

  // compute derivatives for the second part
  double value2 = 0;
  LocalDerivativeType derivative2;
  derivative2.Fill(NumericTraits<LocalDerivativeValueType>::ZeroValue());

  for (MovingPointIterator it = this->m_TransformedPointSet->GetPoints()->Begin(); it != this->m_TransformedPointSet->GetPoints()->End(); ++it) {
    const typename MovingPointSetType::PointType transformedPoint = it.Value();
    const double distance = point.SquaredEuclideanDistanceTo(transformedPoint);
    const double expval = std::exp(-distance / scale2);
    value2 += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) {
      derivative2[dim] += expval * (point[dim] - transformedPoint[dim]);
    }
  }

  // local value
  value = value2 / factor2 - 2.0 * value1 / factor1;

  // local derivatives
  for (size_t dim = 0; dim < this->PointDimension; ++dim) {
    derivative[dim] = -2.0 * (derivative2[dim] / (scale2 * factor2) - 2.0 * derivative1[dim] / (scale1 * factor1));
  }
}
}

#endif
