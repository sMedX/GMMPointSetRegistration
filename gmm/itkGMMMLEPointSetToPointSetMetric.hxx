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
  this->m_UseFixedPointSetKdTree  = true;
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMMLEPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();

  this->InitializeTransformedMovingPointSet();

  this->m_NormalizingValueFactor = -2.0 / (this->m_MovingPointSet->GetNumberOfPoints() * this->m_FixedPointSet->GetNumberOfPoints());

  this->m_NormalizingDerivativeFactor = 4.0 / (this->m_MovingPointSet->GetNumberOfPoints() * this->m_FixedPointSet->GetNumberOfPoints() * this->m_Scale * this->m_Scale);
}

/** Initialize data for current iteration with the input parameters */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMMLEPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::InitializeForIteration(const ParametersType & parameters) const
{
  this->m_Transform->SetParameters(parameters);

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) 
  {
    m_TransformedMovingPointSet->GetPoints()->SetElement(it.Index(), m_Transform->TransformPoint(it.Value()));
  }
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointIterator & it) const
{
  const MovingPointType point = this->m_TransformedMovingPointSet->GetPoints()->at(it.Index());

  const double scale = this->m_Scale * this->m_Scale;

  MeasureType value = NumericTraits<MeasureType>::ZeroValue();

  FixedNeighborsIdentifierType idx;
  this->m_FixedPointsLocator->Search(point, this->m_SearchRadius * this->m_Scale, idx);

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) 
  {
    const FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / scale);
    value += expval;
  }

  return value;
}

template<typename TFixedPointSet, typename TMovingPointSet>
bool
GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative(const MovingPointIterator & it, MeasureType & value, LocalDerivativeType & derivative) const
{
  const MovingPointType point = this->m_TransformedMovingPointSet->GetPoints()->at(it.Index());

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
    const FixedPointType fixedPoint = this->GetFixedPoint(*it);

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
