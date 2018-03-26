#ifndef itkGMMMLEPointSetToPointSetMetric_hxx
#define itkGMMMLEPointSetToPointSetMetric_hxx

#include "itkGMMMLEPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMMMLEPointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMMLEPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();
  this->SetUseFixedPointSetKdTree(true);
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointType & point) const
{
  const double scale = this->m_Scale * this->m_Scale;

  MeasureType value = NumericTraits<MeasureType>::ZeroValue();

  FixedNeighborsIdentifierType idx;
  this->m_FixedPointsLocator->Search(point, this->m_SearchRadius * this->m_Scale, idx);

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) {
    const double distance = point.SquaredEuclideanDistanceTo(this->m_FixedPointSet->GetPoint(*it));
    const double expval = std::exp(-distance / scale);
    value += expval;
  }

  return value;
}

template<typename TFixedPointSet, typename TMovingPointSet>
void
GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative(const MovingPointType & point, MeasureType & value, LocalDerivativeType & derivative) const
{
  const double scale = this->m_Scale * this->m_Scale;

  value = NumericTraits<MeasureType>::ZeroValue();

  derivative.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  if (this->m_UseFixedPointSetKdTree) {
    FixedNeighborsIdentifierType idx;
    this->m_FixedPointsLocator->Search(point, this->m_SearchRadius * this->m_Scale, idx);

    for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) {
      const FixedPointType & fixedPoint = this->m_FixedPointSet->GetPoint(*it);
      const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
      const double expval = std::exp(-distance / scale);
      value += expval;

      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        derivative[dim] += expval * (point[dim] - fixedPoint[dim]);
      }
    }
  }
  else {
    for (FixedPointIterator it = this->m_FixedPointSet->GetPoints()->Begin(); it != this->m_FixedPointSet->GetPoints()->End(); ++it) {
      const FixedPointType & fixedPoint = it.Value();
      const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
      const double expval = std::exp(-distance / scale);
      value += expval;

      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        derivative[dim] += expval * (point[dim] - fixedPoint[dim]);
      }
    }
  }
}
}

#endif
