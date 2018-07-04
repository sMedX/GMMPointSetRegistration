#ifndef itkGMMScalePointSetMetric_h
#define itkGMMScalePointSetMetric_h

#include "itkSingleValuedCostFunction.h"

namespace itk
{
/** \class GMMSigmaMtric
 * \brief Computes similarity between pixel values of a point set and
 * intensity values of an image.
 *
 * This metric computes the average squared differences between pixels
 * in the point set and the transformed set of pixels.
 *
 * Spatial correspondence between both images is established through a
 * Transform.
 */
template< typename TFixedPointSet, typename TMovingPointSet = TFixedPointSet >
class GMMScalePointSetMetric : public SingleValuedCostFunction
{
public:
  /** Standard class typedefs. */
  typedef GMMScalePointSetMetric              Self;
  typedef SingleValuedCostFunction    Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GMMScalePointSetMetric, SingleValuedCostFunction);

  /**  Type of the fixed point set. */
  typedef GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet>  PointSetMetricType;
  typedef typename Superclass::MeasureType                                  MeasureType;
  typedef typename Superclass::ParametersType                               ParametersType;
  typedef typename Superclass::DerivativeType                               DerivativeType;
  typedef typename DerivativeType::ValueType                                DerivativeValueType;

  /** Get the value for single valued optimizers. */
  MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const ParametersType & parameters, DerivativeType & Derivative) const ITK_OVERRIDE;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const ParametersType & parameters, MeasureType & Value, DerivativeType & Derivative) const ITK_OVERRIDE;

  /** Initialize the Metric by making sure that all the components are present and plugged together correctly     */
  virtual void Initialize(void)  throw (ExceptionObject);

  itkSetObjectMacro(PointSetMetric, PointSetMetricType);
  itkGetObjectMacro(PointSetMetric, PointSetMetricType);

  /** Return the number of parameters required by the Transform */
  virtual unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE
  {
    return m_NumberOfParameters;
  }

protected:
  GMMScalePointSetMetric();
  virtual ~GMMScalePointSetMetric() {}

  mutable typename PointSetMetricType::Pointer m_PointSetMetric;
  size_t m_NumberOfParameters;

private:
  GMMScalePointSetMetric(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGMMScalePointSetMetric.hxx"
#endif

#endif
