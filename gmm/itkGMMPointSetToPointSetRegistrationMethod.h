#ifndef itkGMMPointSetToPointSetRegistrationMethod_h
#define itkGMMPointSetToPointSetRegistrationMethod_h

#include "itkProcessObject.h"
#include "itkPointSetToPointSetMetric.h"
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkDataObjectDecorator.h"
#include "itkGMMPointSetToPointSetMetricBase.h"
#include "itkGMMScalePointSetMetric.h"
#include "itkGMMScalePointSetMetricEstimator.h"

namespace itk
{
template< typename TFixedPointSet, typename TMovingPointSet = TFixedPointSet >
class GMMPointSetToPointSetRegistrationMethod : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef GMMPointSetToPointSetRegistrationMethod  Self;
  typedef ProcessObject                                Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GMMPointSetToPointSetRegistrationMethod, ProcessObject);

  /**  Type of the Fixed PointSet. */
  typedef          TFixedPointSet                           FixedPointSetType;
  typedef typename FixedPointSetType::Pointer               FixedPointSetPointer;
  typedef typename FixedPointSetType::ConstPointer          FixedPointSetConstPointer;
  typedef typename FixedPointSetType::PointsContainer       FixedPointsContainerType;
  typedef typename FixedPointsContainerType::ConstIterator  FixedPointConstIterator;

  /**  Type of the Moving PointSet. */
  typedef          TMovingPointSet                          MovingPointSetType;
  typedef typename MovingPointSetType::Pointer              MovingPointSetPointer;
  typedef typename MovingPointSetType::ConstPointer         MovingPointSetConstPointer;
  typedef typename MovingPointSetType::PointsContainer      MovingPointsContainerType;
  typedef typename MovingPointsContainerType::ConstIterator MovingPointConstIterator;

  /**  Type of the metric. */
  typedef GMMPointSetToPointSetMetricBase<FixedPointSetType, MovingPointSetType>  MetricType;
  typedef typename MetricType::Pointer                                            MetricPointer;
  typedef itk::Array<typename MetricType::MeasureType>                            MetricValuesType;
  typedef itk::Array<double>                                                      ScaleType;
  typedef itk::GMMScalePointSetMetricEstimator<MetricType>                        MetricEstimatorType;

  /**  Type of the Transform . */
  typedef typename MetricType::TransformType TransformType;
  typedef typename TransformType::Pointer    TransformPointer;

  /** Type for the output: Using Decorator pattern for enabling
   *  the Transform to be passed in the data pipeline */
  typedef  DataObjectDecorator< TransformType >      TransformOutputType;
  typedef typename TransformOutputType::Pointer      TransformOutputPointer;
  typedef typename TransformOutputType::ConstPointer TransformOutputConstPointer;

  /**  Type of the optimizer. */
  typedef SingleValuedNonLinearOptimizer OptimizerType;
  typedef typename OptimizerType::Pointer OptimizerPointer;

  /** Type of the Transformation parameters This is the same type used to
   *  represent the search space of the optimization algorithm */
  typedef typename MetricType::TransformParametersType  ParametersType;
  typedef typename ParametersType::ValueType            ParametersValueType;

  /** Smart Pointer type to a DataObject. */
  typedef typename DataObject::Pointer DataObjectPointer;

  /** Set/Get the Fixed PointSet. */
  itkSetConstObjectMacro(FixedPointSet, FixedPointSetType);
  itkGetConstObjectMacro(FixedPointSet, FixedPointSetType);

  /** Set/Get the Moving PointSet. */
  itkSetConstObjectMacro(MovingPointSet, MovingPointSetType);
  itkGetConstObjectMacro(MovingPointSet, MovingPointSetType);

  /** Set/Get the Optimizer. */
  itkSetObjectMacro(Optimizer,  OptimizerType);
  itkGetModifiableObjectMacro(Optimizer, OptimizerType);

  /** Set/Get the Metric. */
  itkSetObjectMacro(Metric, MetricType);
  itkGetModifiableObjectMacro(Metric, MetricType);

  /** Set/Get the Metric Estimator. */
  itkSetObjectMacro(MetricEstimator, MetricEstimatorType);
  itkGetModifiableObjectMacro(MetricEstimator, MetricEstimatorType);

  /** Set/Get the Transform. */
  itkSetObjectMacro(Transform, TransformType);
  itkGetModifiableObjectMacro(Transform, TransformType);

  /** Set/Get the initial transformation parameters. */
  virtual void SetInitialTransformParameters(const ParametersType & param);

  itkGetConstReferenceMacro(InitialTransformParameters, ParametersType);

  /** Get the last transformation parameters visited by the optimizer. */
  itkGetConstReferenceMacro(FinalTransformParameters, ParametersType);

  /** Initialize by setting the interconnects between the components. */
  void Initialize() throw (ExceptionObject);

  /** Returns the transform resulting from the registration process  */
  const TransformOutputType * GetOutput() const;

  /** Make a DataObject of the correct type to be used as the specified output. */
  typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

  /** Method to return the latest modified time of this object or any of its cached ivars */
  virtual ModifiedTimeType GetMTime() const ITK_OVERRIDE;

  itkSetMacro(NumberOfLevels, size_t);
  itkGetMacro(NumberOfLevels, size_t);

  /** Set/Get the gradient convergence tolerance. This is a positive real number that determines the accuracy with which the solution is to
  * be found. The optimization terminates when: ||G|| < gtol max(1,||X||) where ||.|| denotes the Euclidean norm.  */
  itkSetMacro(GradientConvergenceTolerance, double);
  itkGetMacro(GradientConvergenceTolerance, double);

  itkGetMacro(InitialMetricValues, MetricValuesType);
  itkGetMacro(FinalMetricValues, MetricValuesType);
  itkGetMacro(MetricParameters, ScaleType);

protected:
  GMMPointSetToPointSetRegistrationMethod();
  virtual ~GMMPointSetToPointSetRegistrationMethod() {};

  /** Method invoked by the pipeline in order to trigger the computation of the registration. */
  virtual void GenerateData() ITK_OVERRIDE;
  double ComputeScaleForPointSetMetric();

  FixedPointSetConstPointer m_FixedPointSet;
  MovingPointSetConstPointer m_MovingPointSet;
  MetricPointer m_Metric;
  OptimizerPointer m_Optimizer;
  typename MetricEstimatorType::Pointer m_MetricEstimator;

  TransformPointer m_Transform;
  ParametersType m_InitialTransformParameters;
  ParametersType m_FinalTransformParameters;

  MetricValuesType m_InitialMetricValues;
  MetricValuesType m_FinalMetricValues;

  ScaleType m_MetricParameters;

  size_t m_NumberOfLevels;
  double m_GradientConvergenceTolerance;

private:
  GMMPointSetToPointSetRegistrationMethod(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGMMPointSetToPointSetRegistrationMethod.hxx"
#endif

#endif
