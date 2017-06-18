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
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMML2RigidPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();
}

/**
 * Get the match Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{
  itkExceptionMacro(<< "not implemented");
}

/**
 * Get the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/*
 * Get both the match Measure and the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  this->m_Transform->SetParameters(parameters);

  if (derivative.size() != this->m_NumberOfParameters) {
    derivative.set_size(this->m_NumberOfParameters);
  }

  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());
  value = NumericTraits<MeasureType>::ZeroValue();

  GradientType gradient;

  double scale = 0.5 * (this->m_FixedPointSetScale*this->m_FixedPointSetScale + this->m_MovingPointSetScale*this->m_MovingPointSetScale);

  for (MovingPointIterator movingIter = this->m_MovingPointSet->GetPoints()->Begin(); movingIter != this->m_MovingPointSet->GetPoints()->End(); ++movingIter) {
    const typename MovingPointSetType::PointType transformedPoint = this->m_Transform->TransformPoint(movingIter.Value());

    gradient.Fill(0);

    for (FixedPointIterator fixedIter = this->m_FixedPointSet->GetPoints()->Begin(); fixedIter != this->m_FixedPointSet->GetPoints()->End(); ++fixedIter) {
      const typename FixedPointSetType::PointType fixedPoint = fixedIter.Value();
      const double distance = transformedPoint.SquaredEuclideanDistanceTo(fixedPoint);
      const double expval = exp(-distance / scale);
      value += expval;

      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        gradient[dim] += (-2.0) * expval * (transformedPoint[dim] - fixedPoint[dim]);
      }
    }

    // compute the derivatives
    this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(movingIter.Value(), this->m_Jacobian, this->m_JacobianCache);

    for (size_t par = 0; par < this->m_NumberOfParameters; par++) {
      for (size_t dim = 0; dim < Self::FixedPointSetDimension; dim++) {
        derivative[par] += this->m_Jacobian(dim, par) * gradient[dim];
      }
    }
  }

  const double factor = this->m_MovingPointSet->GetNumberOfPoints() * this->m_FixedPointSet->GetNumberOfPoints();

  value *= -2.0 / factor;

  for (size_t par = 0; par < this->m_NumberOfParameters; par++) {
    derivative[par] = (-2.0) * derivative[par] / (scale * factor);
  }
}
}

#endif
