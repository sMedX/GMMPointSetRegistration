#ifndef itkGMMKCKdTreePointSetToPointSetMetric_hxx
#define itkGMMKCPointSetToPointSetMetric_hxx

#include "itkGMMKCKdTreePointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMKCKdTreePointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMMKCKdTreePointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMKCKdTreePointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize(void)
throw (ExceptionObject)
{
  Superclass::Initialize();

  typename FixedPointSetType::PointsContainer::ConstPointer fixedPointContainer = m_FixedPointSet->GetPoints();

  m_FixedPointsLocator = FixedPointsLocatorType::New();
  m_FixedPointsLocator->SetPoints(const_cast<typename FixedPointSetType::PointsContainer*> (fixedPointContainer.GetPointer()));
  m_FixedPointsLocator->Initialize();
}

/**
 * Get the match Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMMKCKdTreePointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMMKCKdTreePointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{
  itkExceptionMacro(<< "not implemented");
}
/**
 * Get the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMKCKdTreePointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMKCKdTreePointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  m_Transform->SetParameters(parameters);

  if (derivative.size() != m_NumberOfParameters) {
    derivative.set_size(m_NumberOfParameters);
  }

  m_TransformedPointSet = MovingPointSetType::New();
  for (MovingPointIterator iter = m_MovingPointSet->GetPoints()->Begin(); iter != m_MovingPointSet->GetPoints()->End(); ++iter) {
    m_TransformedPointSet->SetPoint(iter.Index(), m_Transform->TransformPoint(iter.Value()));
  }

  typename MovingPointsLocatorType::Pointer m_TransformedPointsLocator = MovingPointsLocatorType::New();
  m_TransformedPointsLocator->SetPoints(m_TransformedPointSet->GetPoints());
  m_TransformedPointsLocator->Initialize();

  double value1 = 0;
  double value2 = 0;

  GradientType gradient1;
  GradientType gradient2;

  LocalDerivativeType derivative1(m_NumberOfParameters);
  derivative1.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  LocalDerivativeType derivative2(m_NumberOfParameters);
  derivative2.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  const double scale1 = (m_FixedPointSetScale*m_FixedPointSetScale + m_MovingPointSetScale*m_MovingPointSetScale) / 2;
  const double scale2 = m_MovingPointSetScale*m_MovingPointSetScale;

  const double radius1 = m_Radius * sqrt(scale1);
  const double radius2 = m_Radius * sqrt(scale2);

  for (MovingPointIterator movingIter1 = m_TransformedPointSet->GetPoints()->Begin(); movingIter1 != m_TransformedPointSet->GetPoints()->End(); ++movingIter1) {
    const typename MovingPointSetType::PointType transformedPoint1 = movingIter1.Value();

    //------------------------------------
    // compute gradient for the first part
    gradient1.Fill(0);

    typename FixedPointsLocatorType::NeighborsIdentifierType fixedPointNeighbors;
    m_FixedPointsLocator->FindPointsWithinRadius(transformedPoint1, radius1, fixedPointNeighbors);

    for (size_t n = 0; n < fixedPointNeighbors.size(); ++n) {
      const typename FixedPointSetType::PointType fixedPoint = m_FixedPointSet->GetPoint(fixedPointNeighbors[n]);
      const double distance = transformedPoint1.SquaredEuclideanDistanceTo(fixedPoint) / scale1;
      value1 += exp(-distance);

      for (size_t dim = 0; dim < PointDimension; ++dim) {
        gradient1[dim] += (-2.0) * exp(-distance) * (transformedPoint1[dim] - fixedPoint[dim]) / scale1;
      }
    }

    //-------------------------------------
    // compute gradient for the second part
    gradient2.Fill(0);

    typename MovingPointsLocatorType::NeighborsIdentifierType transformedPointNeighbors;
    m_TransformedPointsLocator->FindPointsWithinRadius(transformedPoint1, radius2, transformedPointNeighbors);

    for (size_t n = 0; n < transformedPointNeighbors.size(); ++n) {
      const typename MovingPointSetType::PointType transformedPoint2 = m_TransformedPointSet->GetPoint(transformedPointNeighbors[n]);
      const double distance = transformedPoint1.SquaredEuclideanDistanceTo(transformedPoint2) / scale2;
      value2 += exp(-distance);

      for (size_t dim = 0; dim < PointDimension; ++dim) {
        gradient2[dim] += (-2.0) * exp(-distance) * (transformedPoint1[dim] - transformedPoint2[dim]) / scale2;
      }
    }

    //-------------------------------------
    // compute derivatives
    m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(m_MovingPointSet->GetPoint(movingIter1.Index()), m_Jacobian, m_JacobianCache);

    for (size_t par = 0; par < m_NumberOfParameters; par++) {
      for (size_t dim = 0; dim < PointDimension; dim++) {
        derivative1[par] += m_Jacobian(dim, par) * gradient1[dim];
        derivative2[par] += m_Jacobian(dim, par) * gradient2[dim];
      }
    }
  }

  const double factor = m_MovingPointSet->GetNumberOfPoints() * m_FixedPointSet->GetNumberOfPoints();
  const double ratio = value1 / value2;
  
  value = -value1 * ratio / factor;

  for (size_t par = 0; par < m_NumberOfParameters; par++) {
    derivative[par] = (-2.0) * (derivative1[par] - derivative2[par] * ratio) * ratio / factor;
  }
}
}

#endif
