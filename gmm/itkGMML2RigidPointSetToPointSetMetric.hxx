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

  this->m_NormalizingValueFactor = -2.0 / (this->m_MovingPointSet->GetNumberOfPoints() * this->m_FixedPointSet->GetNumberOfPoints());

  this->m_NormalizingDerivativeFactor = -2.0 * this->m_NormalizingValueFactor / (this->m_Scale * this->m_Scale);
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointType & point) const
{
  const double scale = this->m_Scale * this->m_Scale;

  MeasureType value = NumericTraits<MeasureType>::ZeroValue();

  if (this->m_UseFixedPointSetKdTree) {
    FixedNeighborsIdentifierType idx;
    this->m_FixedPointsLocator->Search(point, this->m_SearchRadius * this->m_Scale, idx);

    for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) {
      const double distance = point.SquaredEuclideanDistanceTo(this->m_FixedPointSet->GetPoint(*it));
      const double expval = std::exp(-distance / scale);
      value += expval;
    }
  }
  else {
    for (FixedPointIterator it = this->m_FixedPointSet->GetPoints()->Begin(); it != this->m_FixedPointSet->GetPoints()->End(); ++it) {
      const double distance = point.SquaredEuclideanDistanceTo(it.Value());
      const double expval = std::exp(-distance / scale);
      value += expval;
    }
  }

  return value;
}

template<typename TFixedPointSet, typename TMovingPointSet>
void
GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
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
