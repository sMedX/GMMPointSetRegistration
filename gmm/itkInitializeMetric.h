#pragma once

#include "itkGMML2RigidPointSetToPointSetMetric.h"
#include "itkGMML2PointSetToPointSetMetric.h"
#include "itkGMMKCPointSetToPointSetMetric.h"
#include "itkGMMMLEPointSetToPointSetMetric.h"

namespace itk
{
  template <typename TFixedPointSet, typename TMovingPointSet>
  class InitializeMetric : public itk::ProcessObject
  {
  public:
    /** Standard class typedefs. */
    typedef InitializeMetric                        Self;
    typedef itk::ProcessObject                      Superclass;
    typedef itk::SmartPointer<Self>                 Pointer;
    typedef itk::SmartPointer<const Self>           ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    itkTypeMacro(InitializeMetric, itk::ProcessObject);

    enum class Metric
    {
      L2Rigid,
      L2,
      KC,
      MLE,
    };

    /** typedefs */
    typedef itk::GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet> MetricType;

    // Get metric
    itkGetObjectMacro(Metric, MetricType);

    // Set/Get type of metric
    itkSetEnumMacro(TypeOfMetric, Metric);
    itkGetEnumMacro(TypeOfMetric, Metric);

    void SetTypeOfMetric(const size_t & type) 
    { 
      this->SetTypeOfMetric(static_cast<Metric>(type)); 
    }

    void Update()
    {
      switch (m_TypeOfMetric) {
      case Metric::L2Rigid: {
        typedef itk::GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> GMML2RigidMetricType;
        m_Metric = GMML2RigidMetricType::New();
        break;
      }
      case Metric::L2:{
        typedef itk::GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> GMML2MetricType;
        m_Metric = GMML2MetricType::New();
        break;
      }
      case Metric::KC: {
        typedef itk::GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> GMMKCMetricType;
        m_Metric = GMMKCMetricType::New();
        break;
      }
      case Metric::MLE:{
        typedef itk::GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> GMMMLEMetricType;
        m_Metric = GMMMLEMetricType::New();
        break;
      }
      default:
        itkExceptionMacro(<< "Invalid type of metric");
      }
    }

    void PrintReport() const
    {
      std::cout << "class name " << this->GetNameOfClass() << std::endl;
      std::cout << "metric     " << m_Metric->GetNameOfClass() << std::endl;
      std::cout << std::endl;
    }

  protected:
    Metric m_TypeOfMetric = Metric::L2Rigid;
    typename MetricType::Pointer m_Metric = ITK_NULLPTR;

    InitializeMetric()
    {
      this->SetNumberOfRequiredInputs(0);
      this->SetNumberOfRequiredOutputs(0);
    }
    ~InitializeMetric() {}
  };
}
