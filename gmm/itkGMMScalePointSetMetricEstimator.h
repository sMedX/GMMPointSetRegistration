#ifndef itkGMMScalePointSetMetricEstimator_h
#define itkGMMScalePointSetMetricEstimator_h

#include "itkObject.h"
#include "itkGMMScalePointSetMetric.h"

namespace itk
{
template< typename TMetricType >
class GMMScalePointSetMetricEstimator : public Object
{
public:
  /** Standard class typedefs. */
  typedef GMMScalePointSetMetricEstimator  Self;
  typedef Object                                       Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GMMScalePointSetMetricEstimator, Superclass);

  typedef TMetricType MetricType;
  typedef typename MetricType::ParametersType                            ParametersType;

  typedef typename MetricType::FixedPointSetType                         FixedPointSetType;
  typedef typename FixedPointSetType::PointsContainer::ConstIterator     FixedPointIterator;

  typedef typename MetricType::MovingPointSetType                        MovingPointSetType;
  typedef typename MovingPointSetType::PointsContainer::ConstIterator    MovingPointIterator;

  /** Method invoked by the pipeline in order to trigger the computation of the registration. */
  virtual void Estimate() ITK_OVERRIDE;

  itkSetObjectMacro(Metric, MetricType);
  itkGetObjectMacro(Metric, MetricType);

  virtual void SetInitialParameters(const ParametersType &);
  itkGetMacro(InitialParameters, ParametersType);
  itkGetMacro(Parameters, ParametersType);

  itkSetMacro(UseInitialParameters, bool);
  itkGetMacro(UseInitialParameters, bool);

  itkSetMacro(Trace, bool);
  itkGetMacro(Trace, bool);

protected:
  GMMScalePointSetMetricEstimator();
  virtual ~GMMScalePointSetMetricEstimator() {};

  void ComputeDistance();

  typename MetricType::Pointer m_Metric;
  ParametersType m_InitialParameters;
  ParametersType m_Parameters;
  bool m_UseInitialParameters;

  bool m_Trace;
  double m_RMSDistance;

private:
  GMMScalePointSetMetricEstimator(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGMMScalePointSetMetricEstimator.hxx"
#endif

#endif
