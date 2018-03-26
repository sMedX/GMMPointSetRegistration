#pragma once

#include "itkGMML2RigidPointSetToPointSetMetric.h"
#include "itkGMML2PointSetToPointSetMetric.h"
#include "itkGMMKCPointSetToPointSetMetric.h"
#include "itkGMMMLEPointSetToPointSetMetric.h"

namespace itk
{
  template <typename TFixedPointSet, typename TMovingPointSet>
  class InitializeMetric : public itk::Object
  {
  public:
    /** Standard class typedefs. */
    typedef InitializeMetric                        Self;
    typedef itk::Object                             Superclass;
    typedef itk::SmartPointer<Self>                 Pointer;
    typedef itk::SmartPointer<const Self>           ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    itkTypeMacro(InitializeMetric, Object);

    /** typedefs */
    typedef itk::GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet> MetricType;

    // define type of metric
    enum class Metric
    {
      GMML2Rigid = 0,
      GMML2 = 1,
      GMMKC = 2,
      GMMMLE = 3
    };

    itkSetEnumMacro(TypeOfMetric, Metric);
    itkGetEnumMacro(TypeOfMetric, Metric);
    void SetTypeOfMetric(const size_t & type) { this->SetTypeOfMetric(static_cast<Metric>(type)); }

    // Get metric
    itkGetObjectMacro(Metric, MetricType);

    void Initialize()
    {
      switch (m_TypeOfMetric) {
      case Metric::GMML2Rigid: {
        typedef itk::GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> MetricType;
        m_Metric = MetricType::New();
        break;
      }
      case Metric::GMML2:{
        typedef itk::GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> MetricType;
        m_Metric = MetricType::New();
        break;
      }
      case Metric::GMMKC: {
        typedef itk::GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> MetricType;
        m_Metric = MetricType::New();
        break;
      }
      case Metric::GMMMLE: {
        typedef itk::GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> MetricType;
        m_Metric = MetricType::New();
        break;
      }
      default: {
        itkExceptionMacro(<< "Unknown type of the metric to initialize");
        return;
      }
      }

      if (m_Metric == nullptr) {
        itkExceptionMacro(<< "metric has not been initialized.");
      }
    }

    void PrintReport() const
    {
      std::cout << "class name " << this->GetNameOfClass() << std::endl;
      std::cout << "metric     " << m_Metric->GetNameOfClass() << std::endl;
      std::cout << std::endl;
    }

  protected:
    Metric m_TypeOfMetric;
    typename MetricType::Pointer m_Metric;

    InitializeMetric() 
    {
      m_Metric = nullptr;
    }
    ~InitializeMetric() {}
  };
}
