#ifndef itkGMML2RigidPointSetToPointSetMetric_hxx
#define itkGMML2RigidPointSetToPointSetMetric_hxx

#include "itkGMML2RigidPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMML2RigidPointSetToPointSetMetric()
{
  this->m_UseFixedPointSetKdTree = true;
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMML2RigidPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();

  this->m_NormalizingValueFactor = -1.0 / (this->m_FixedPointSet->GetNumberOfPoints() * this->m_MovingPointSet->GetNumberOfPoints());
  this->m_NormalizingDerivativeFactor = -2.0 * this->m_NormalizingValueFactor / this->m_Variance;
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointIterator & it) const
{
  const FixedPointType point = this->m_Transform->TransformPoint(it.Value());

  MeasureType value = NumericTraits<MeasureType>::ZeroValue();

  FixedNeighborsIdentifierType idx;
  this->SearchFixedPoints(point, idx);

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) 
  {
    FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / this->m_Variance);
    value += expval;
  }

  return value;
}

template<typename TFixedPointSet, typename TMovingPointSet>
bool
GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative(const MovingPointIterator & it, MeasureType & value, LocalDerivativeType & derivative) const
{
  const FixedPointType point = this->m_Transform->TransformPoint(it.Value());

  FixedNeighborsIdentifierType idx;
  if (!this->SearchFixedPoints(point, idx)) 
  {
    return false;
  }

  value = NumericTraits<MeasureType>::ZeroValue();
  derivative.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) 
  {
    FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / this->m_Variance);
    value += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) 
    {
      derivative[dim] += expval * (point[dim] - fixedPoint[dim]);
    }
  }

  return true;
}

template<typename TFixedPointSet, typename TMovingPointSet>
bool
GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodDerivative(const MovingPointIterator & it, LocalDerivativeType & derivative) const
{
  const FixedPointType point = this->m_Transform->TransformPoint(it.Value());

  FixedNeighborsIdentifierType idx;
  if (!this->SearchFixedPoints(point, idx)) {
    return false;
  }

  derivative.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) {
    FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / this->m_Variance);

    for (size_t dim = 0; dim < this->PointDimension; ++dim) {
      derivative[dim] += expval * (point[dim] - fixedPoint[dim]);
    }
  }

  return true;
}
}

#endif
