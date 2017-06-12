#pragma once

#include "itkGMMPointSetToPointSetMetricBase.h"
#include "gmm/itkGMMRigidPointSetToPointSetMetric.h"
#include "gmm/itkGMMPointSetToPointSetMetric.h"
#include "gmm/itkGMMKCPointSetToPointSetMetric.h"
#include "gmm/itkGMMKCKdTreePointSetToPointSetMetric.h"
#include "gmm/itkGMMMLEPointSetToPointSetMetric.h"

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
      Rigid,
      GMM,
      KC,
      KCKdTree,
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
      case Metric::Rigid: {
        typedef itk::GMMRigidPointSetToPointSetMetric<PointSetType, PointSetType> MetricType;
        m_Metric = MetricType::New();
        break;
      }
      case Metric::GMM:{
        typedef itk::GMMPointSetToPointSetMetric<PointSetType, PointSetType> MetricType;
        m_Metric = MetricType::New();
        break;
      }
      case Metric::KC: {
        typedef itk::GMMKCPointSetToPointSetMetric<PointSetType, PointSetType> MetricType;
        m_Metric = MetricType::New();
        break;
      }
      case Metric::KCKdTree: {
        typedef itk::GMMKCKdTreePointSetToPointSetMetric<PointSetType, PointSetType> MetricType;
        m_Metric = MetricType::New();
        break;
      }
      case Metric::MLE:{
        typedef itk::GMMMLEPointSetToPointSetMetric<PointSetType, PointSetType> MetricType;
        m_Metric = MetricType::New();
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
    Metric m_TypeOfMetric = Metric::Rigid;
    typename MetricType::Pointer m_Metric = nullptr;

    InitializeMetric()
    {
      this->SetNumberOfRequiredInputs(0);
      this->SetNumberOfRequiredOutputs(0);
    }
    ~InitializeMetric() {}
  };
}
