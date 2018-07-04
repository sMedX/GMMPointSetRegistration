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
  this->m_UseFixedPointSetKdTree = true;
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMML2PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();

  this->InitializeTransformedMovingPointSet();

  this->m_NormalizingValueFactor = 1.0 / this->m_TransformedMovingPointSet->GetNumberOfPoints();

  this->m_NormalizingDerivativeFactor = -2.0 / (this->m_TransformedMovingPointSet->GetNumberOfPoints() * this->m_Scale * this->m_Scale);
}

/** Initialize data for current iteration with the input parameters */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMML2PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::InitializeForIteration(const ParametersType & parameters) const
{
  this->m_Transform->SetParameters(parameters);

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) 
  {
    m_TransformedMovingPointSet->GetPoints()->SetElement(it.Index(), m_Transform->TransformPoint(it.Value()));
  }
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointIterator & it) const
{
  const MovingPointType point = this->GetTransformedMovingPoint(it.Index());

  const double scale = this->m_Scale * this->m_Scale;

  // compute value for the first sum
  double value1 = 0;

  FixedNeighborsIdentifierType idx;
  this->SearchFixedPoints(point, idx);

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) 
  {
    const FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / scale);
    value1 += expval;
  }

  // compute value for the second sum
  double value2 = 0;
 
  for (FixedPointIterator it = this->m_TransformedMovingPointSet->GetPoints()->Begin(); it != this->m_TransformedMovingPointSet->GetPoints()->End(); ++it) 
  {
    const double distance = point.SquaredEuclideanDistanceTo(it.Value());
    const double expval = std::exp(-distance / scale);
    value2 += expval;
  }

  // local value
  const double factor1 = this->m_FixedPointSet->GetNumberOfPoints();
  const double factor2 = this->m_MovingPointSet->GetNumberOfPoints();

  const double value = -2.0 * value1 / factor1 + value2 / factor2;

  return value;
}

template<typename TFixedPointSet, typename TMovingPointSet>
bool
GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative(const MovingPointIterator & it, MeasureType & value, LocalDerivativeType & derivative) const
{
  const MovingPointType point = this->GetTransformedMovingPoint(it.Index());

  const double scale = this->m_Scale * this->m_Scale;

  // compute value and derivative gradient for the first sum
  double value1 = 0;
  LocalDerivativeType derivative1;
  derivative1.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  FixedNeighborsIdentifierType idx;
  this->SearchFixedPoints(point, idx);

  for (FixedNeighborsIteratorType it = idx.begin(); it != idx.end(); ++it) 
  {
    const FixedPointType fixedPoint = this->GetFixedPoint(*it);

    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / scale);
    value1 += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) 
    {
      derivative1[dim] += expval * (point[dim] - fixedPoint[dim]);
    }
  }

  // compute derivatives for the second part
  double value2 = 0;
  LocalDerivativeType derivative2;
  derivative2.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  for (FixedPointIterator it = this->m_TransformedMovingPointSet->GetPoints()->Begin(); it != this->m_TransformedMovingPointSet->GetPoints()->End(); ++it) 
  {
    const FixedPointType transformedPoint = it.Value();
    const double distance = point.SquaredEuclideanDistanceTo(transformedPoint);
    const double expval = std::exp(-distance / scale);
    value2 += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) 
    {
      derivative2[dim] += expval * (point[dim] - transformedPoint[dim]);
    }
  }

  // local value
  const double factor1 = this->m_FixedPointSet->GetNumberOfPoints();
  const double factor2 = this->m_MovingPointSet->GetNumberOfPoints();

  value = -2.0 * value1 / factor1 + value2 / factor2;

  // local derivatives
  for (size_t dim = 0; dim < this->PointDimension; ++dim) 
  {
    derivative[dim] = -2.0 * derivative1[dim] / factor1 + derivative2[dim] / factor2;
  }

  return true;
}
}

#endif
