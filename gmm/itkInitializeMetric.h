#pragma once

#include "itkGMML2RigidPointSetToPointSetMetric.h"
#include "itkGMML2PointSetToPointSetMetric.h"
#include "itkGMMKCPointSetToPointSetMetric.h"

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

    enum class Metric
    {
      GMML2Rigid,
      GMML2,
      GMMKC,
    };

    /** typedefs */
    typedef itk::GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet> MetricType;

    // Get metric
    itkGetObjectMacro(Metric, MetricType);

    // Set/Get type of metric
    itkSetEnumMacro(TypeOfMetric, Metric);
    itkGetEnumMacro(TypeOfMetric, Metric);
    void SetTypeOfMetric(const size_t & type) { this->SetTypeOfMetric(static_cast<Metric>(type)); }

    void Initialize()
    {
      switch (m_TypeOfMetric) {
      case Metric::GMML2Rigid: {
        typedef itk::GMML2RigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> GMML2RigidMetricType;
        m_Metric = GMML2RigidMetricType::New();
        break;
      }
      case Metric::GMML2:{
        typedef itk::GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> GMML2MetricType;
        m_Metric = GMML2MetricType::New();
        break;
      }
      case Metric::GMMKC: {
        typedef itk::GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet> GMMKCMetricType;
        m_Metric = GMMKCMetricType::New();
        break;
      }
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
    typename MetricType::Pointer m_Metric = nullptr;

    InitializeMetric() {}
    ~InitializeMetric() {}
  };
}
