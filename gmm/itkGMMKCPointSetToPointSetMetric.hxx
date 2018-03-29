#ifndef itkGMMKCPointSetToPointSetMetric_hxx
#define itkGMMKCPointSetToPointSetMetric_hxx

#include "itkGMMKCPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMMKCPointSetToPointSetMetric()
{
  this->m_UseFixedPointSetKdTree = true;
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMKCPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();

  this->InitializeTransformedMovingPointSet();

  const double factor = (double) this->m_MovingPointSet->GetNumberOfPoints() / (double) this->m_FixedPointSet->GetNumberOfPoints();

  this->m_NormalizingValueFactor = -factor / (this->m_MovingPointSet->GetNumberOfPoints() * this->m_FixedPointSet->GetNumberOfPoints());

  this->m_NormalizingDerivativeFactor = -4.0 * factor * this->m_NormalizingValueFactor;
}

/** Initialize data for current iteration with the input parameters */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMKCPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::InitializeForIteration(const ParametersType & parameters) const
{
  this->m_Transform->SetParameters(parameters);

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) 
  {
    m_TransformedMovingPointSet->GetPoints()->SetElement(it.Index(), m_Transform->TransformPoint(it.Value()));
  }
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointIterator & it) const
{
  const MovingPointType point = this->m_TransformedMovingPointSet->GetPoints()->at(it.Index());

  const double scale = this->m_Scale*this->m_Scale;

  // compute value for the first sum
  double value1 = 0;

  FixedNeighborsIdentifierType idx;
  this->m_FixedPointsLocator->Search(point, this->m_SearchRadius * this->m_Scale, idx);

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) 
  {
    const FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / scale);
    value1 += expval;
  }

  // compute value for the second sum
  double value2 = 0;

  for (MovingPointIterator it = this->m_TransformedMovingPointSet->GetPoints()->Begin(); it != this->m_TransformedMovingPointSet->GetPoints()->End(); ++it) 
  {
    const double distance = point.SquaredEuclideanDistanceTo(it.Value());
    const double expval = std::exp(-distance / scale);
    value2 += expval;
  }

  // compute local value
  const double ratio = value1 / value2;
  const double value = value1 * ratio;

  return value;
}

template<typename TFixedPointSet, typename TMovingPointSet>
bool
GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative(const MovingPointIterator & it, MeasureType & value, LocalDerivativeType & derivative) const
{
  const MovingPointType point = this->m_TransformedMovingPointSet->GetPoints()->at(it.Index());

  const double scale = this->m_Scale * this->m_Scale;

  // compute gradient for the first sum
  double value1 = 0;
  LocalDerivativeType derivative1;
  derivative1.Fill(0);

  FixedNeighborsIdentifierType idx;
  this->m_FixedPointsLocator->Search(point, this->m_SearchRadius * this->m_Scale, idx);

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) 
  {
    const FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / scale);
    value1 += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) {
      derivative1[dim] += expval * (point[dim] - fixedPoint[dim]) / scale;
    }
  }

  // compute gradient for the second part
  double value2 = 0;
  LocalDerivativeType derivative2;
  derivative2.Fill(0);

  for (MovingPointIterator it = this->m_TransformedMovingPointSet->GetPoints()->Begin(); it != this->m_TransformedMovingPointSet->GetPoints()->End(); ++it) 
  {
    const MovingPointType transformedPoint = it.Value();

    const double distance = point.SquaredEuclideanDistanceTo(transformedPoint);
    const double expval = std::exp(-distance / scale);
    value2 += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) {
      derivative2[dim] += expval * (point[dim] - transformedPoint[dim]) / scale;
    }
  }

  // compute local value
  const double ratio = value1 / value2;
  value = value1 * ratio;

  // compute local derivatives
  for (size_t dim = 0; dim < this->PointDimension; ++dim) 
  {
    derivative[dim] = (derivative1[dim] - derivative2[dim] * ratio) * ratio;
  }

  return true;
}
}

#endif
