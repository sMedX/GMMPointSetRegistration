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

  this->m_NormalizingValueFactor = -2.0 / (this->m_FixedPointSet->GetNumberOfPoints() * this->m_MovingPointSet->GetNumberOfPoints());

  this->m_NormalizingDerivativeFactor = 4.0 / (this->m_FixedPointSet->GetNumberOfPoints() * this->m_MovingPointSet->GetNumberOfPoints() * this->m_Scale * this->m_Scale);
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
  this->m_FixedPointsLocator->Search(point, this->m_SearchRadius * this->m_Scale, idx);

  if (idx.size() == 0) 
  {
    return value;
  }

  const double scale = this->m_Scale * this->m_Scale;

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) 
  {
    FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / scale);
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
  this->m_FixedPointsLocator->Search(point, this->m_SearchRadius * this->m_Scale, idx);

  if (idx.size() == 0) 
  {
    return false;
  }

  const double scale = this->m_Scale * this->m_Scale;

  value = NumericTraits<MeasureType>::ZeroValue();
  derivative.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) 
  {
    FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / scale);
    value += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) 
    {
      derivative[dim] += expval * (point[dim] - fixedPoint[dim]);
    }
  }

  return true;
}
}

#endif
