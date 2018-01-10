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
  this->SetUseFixedPointSetKdTree(true);
  this->SetUseMovingPointSetKdTree(false);
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMML2PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();

  this->m_NormalizingValueFactor = 1.0 / this->m_MovingPointSet->GetNumberOfPoints();

  this->m_NormalizingDerivativeFactor = -2.0 * this->m_NormalizingValueFactor / (this->m_Scale * this->m_Scale);
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointType & point) const
{
  const double factor1 = this->m_TransformedMovingPointSet->GetNumberOfPoints() * this->m_FixedPointSet->GetNumberOfPoints();
  const double factor2 = this->m_TransformedMovingPointSet->GetNumberOfPoints() * this->m_TransformedMovingPointSet->GetNumberOfPoints();
  const double scale = this->m_Scale * this->m_Scale;

  // compute value for the first sum
  double value1 = 0;

  if (this->m_UseFixedPointSetKdTree) {
    FixedNeighborsIdentifierType idx;
    this->m_FixedPointsLocator->Search(point, this->m_Radius * this->m_Scale, idx);

    for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) {
      const double distance = point.SquaredEuclideanDistanceTo(this->m_FixedPointSet->GetPoint(*it));
      const double expval = std::exp(-distance / scale);
      value1 += expval;
    }
  }
  else {
    for (FixedPointIterator it = this->m_FixedPointSet->GetPoints()->Begin(); it != this->m_FixedPointSet->GetPoints()->End(); ++it) {
      const double distance = point.SquaredEuclideanDistanceTo(it.Value());
      const double expval = std::exp(-distance / scale);
      value1 += expval;
    }
  }

  // compute value for the second sum
  double value2 = 0;
 
  for (MovingPointIterator it = this->m_TransformedMovingPointSet->GetPoints()->Begin(); it != this->m_TransformedMovingPointSet->GetPoints()->End(); ++it) {
    const double distance = point.SquaredEuclideanDistanceTo(it.Value());
    const double expval = std::exp(-distance / scale);
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
  const double factor2 = this->m_TransformedMovingPointSet->GetNumberOfPoints();
  const double scale = this->m_Scale * this->m_Scale;

  // compute value and derivative gradient for the first sum
  double value1 = 0;
  LocalDerivativeType derivative1;
  derivative1.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  if (this->m_UseFixedPointSetKdTree) {
    FixedNeighborsIdentifierType idx;
    this->m_FixedPointsLocator->Search(point, this->m_Radius * this->m_Scale, idx);

    for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) {
      const FixedPointType & fixedPoint = this->m_FixedPointSet->GetPoint(*it);
      const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
      const double expval = std::exp(-distance / scale);
      value1 += expval;

      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        derivative1[dim] += expval * (point[dim] - fixedPoint[dim]);
      }
    }
  }
  else {
    for (FixedPointIterator it = this->m_FixedPointSet->GetPoints()->Begin(); it != this->m_FixedPointSet->GetPoints()->End(); ++it) {
      const FixedPointType & fixedPoint = it.Value();
      const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
      const double expval = std::exp(-distance / scale);
      value1 += expval;

      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        derivative1[dim] += expval * (point[dim] - fixedPoint[dim]);
      }
    }
  }

  // compute derivatives for the second part
  double value2 = 0;
  LocalDerivativeType derivative2;
  derivative2.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  for (MovingPointIterator it = this->m_TransformedMovingPointSet->GetPoints()->Begin(); it != this->m_TransformedMovingPointSet->GetPoints()->End(); ++it) {
    const MovingPointType & transformedPoint = it.Value();
    const double distance = point.SquaredEuclideanDistanceTo(transformedPoint);
    const double expval = std::exp(-distance / scale);
    value2 += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) {
      derivative2[dim] += expval * (point[dim] - transformedPoint[dim]);
    }
  }

  // local value
  value = value2 / factor2 - 2.0 * value1 / factor1;

  // local derivatives
  for (size_t dim = 0; dim < this->PointDimension; ++dim) {
    derivative[dim] = (derivative2[dim] / factor2 - 2.0 * derivative1[dim] / factor1);
  }
}
}

#endif
