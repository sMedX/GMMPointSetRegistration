#pragma once

#include <itkObject.h>
#include <itkListSample.h>
#include <itkSampleToHistogramFilter.h>
#include <itkHistogram.h>
#include <itkPointsLocator.h>

namespace itk
{
  template< typename TFixedPointSet, typename TMovingPointSet = TFixedPointSet >
  class PointSetToPointSetMetrics : public itk::Object
  {
  public:
    /** Standard class typedefs. */
    typedef PointSetToPointSetMetrics               Self;
    typedef itk::Object                             Superclass;
    typedef itk::SmartPointer< Self >               Pointer;
    typedef itk::SmartPointer< const Self >         ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(PointSetToPointSetMetrics, itk::Object);

    /** Types transferred from the base class */
    typedef double                                         MeasureType;
    typedef TFixedPointSet                                 FixedPointSetType;
    typedef TMovingPointSet                                MovingPointSetType;
    typedef typename FixedPointSetType::ConstPointer       FixedPointSetConstPointer;
    typedef typename MovingPointSetType::ConstPointer      MovingPointSetConstPointer;

    /**  Type of the parameters. */
    typedef itk::Statistics::ListSample<typename MovingPointSetType::PointType> ListSampleType;
    typedef typename itk::PointsLocator<typename FixedPointSetType::PointsContainer> FixedPointsLocatorType;
    typedef typename itk::PointsLocator<typename MovingPointSetType::PointsContainer> MovingPointsLocatorType;

    typedef itk::Statistics::ListSample<itk::Vector<MeasureType, 1>> ListMeasureType;
    typedef typename itk::Statistics::Histogram<MeasureType, itk::Statistics::DenseFrequencyContainer2> HistogramType;
    typedef itk::Statistics::SampleToHistogramFilter<ListMeasureType, HistogramType> HistogramFilterType;

    /** Get/Set the Fixed Point Set.  */
    itkSetConstObjectMacro(FixedPointSet, FixedPointSetType);
    itkGetConstObjectMacro(FixedPointSet, FixedPointSetType);

    itkSetConstObjectMacro(TargetPointSet, FixedPointSetType);
    itkGetConstObjectMacro(TargetPointSet, FixedPointSetType);

    /** Get/Set the Moving Image.  */
    itkSetConstObjectMacro(MovingPointSet, MovingPointSetType);
    itkGetConstObjectMacro(MovingPointSet, MovingPointSetType);

    /*Get/Set values to compute quantile. */
    itkSetMacro(LevelOfQuantile, double);
    itkGetMacro(LevelOfQuantile, double);

    itkSetMacro(HistogramSize, size_t);
    itkGetMacro(HistogramSize, size_t);

    /*Get metrics values. */
    itkGetMacro(MeanValue, MeasureType);
    itkGetMacro(RMSEValue, MeasureType);
    itkGetMacro(QuantileValue, MeasureType);
    itkGetMacro(MaximalValue, MeasureType);

    void PrintReport(std::ostream& os) const
    {
      std::string indent = "    ";

      os << "Metric values:" << std::endl;
      os << indent << "    Mean = " << m_MeanValue << std::endl;
      os << indent << "    RMSE = " << m_RMSEValue << std::endl;
      os << indent << "Quantile = " << m_QuantileValue << ", level = " << m_LevelOfQuantile << std::endl;
      os << indent << " Maximal = " << m_MaximalValue << std::endl;
      os << std::endl;

      if (m_TargetPointSet) {
        os << "Target metric values:" << std::endl;
        os << indent << "    Mean = " << m_TargetMeanValue << std::endl;
        os << indent << "    RMSE = " << m_TargetRMSEValue << std::endl;
        os << indent << "Quantile = " << m_TargetQuantileValue << ", level = " << m_LevelOfQuantile << std::endl;
        os << indent << " Maximal = " << m_TargetMaximalValue << std::endl;
        os << std::endl;
      }
    }

    /** Compute metrics. */
    void Compute()
    {
      std::vector<MeasureType> movingToFixedMetrics;
      this->ComputeMetrics<MovingPointSetType, FixedPointSetType>(movingToFixedMetrics, m_MovingPointSet, m_FixedPointSet);

      std::vector<MeasureType> fixedToMovingMetrics;
      this->ComputeMetrics<FixedPointSetType, MovingPointSetType>(fixedToMovingMetrics, m_FixedPointSet, m_MovingPointSet);

      m_MeanValue = 0.5 * (movingToFixedMetrics[0] + fixedToMovingMetrics[0]);
      m_RMSEValue = 0.5 * (movingToFixedMetrics[1] + fixedToMovingMetrics[1]);
      m_QuantileValue = 0.5 * (movingToFixedMetrics[2] + fixedToMovingMetrics[2]);
      m_MaximalValue = 0.5 * (movingToFixedMetrics[3] + fixedToMovingMetrics[3]);

      this->ComputeTargetMetrics();
    }

  protected:
    PointSetToPointSetMetrics() 
    {
      m_FixedPointSet = ITK_NULLPTR;
      m_TargetPointSet = ITK_NULLPTR;
      m_MovingPointSet = ITK_NULLPTR;
    }
    virtual ~PointSetToPointSetMetrics() {}

  private:
    PointSetToPointSetMetrics(const Self &);
    void operator=(const Self &);

    FixedPointSetConstPointer m_FixedPointSet;
    FixedPointSetConstPointer m_TargetPointSet;
    MovingPointSetConstPointer m_MovingPointSet;

    size_t m_BucketSize = 16;
    size_t m_HistogramSize = 1000;
    double m_LevelOfQuantile = 0.95;

    MeasureType m_MeanValue;
    MeasureType m_RMSEValue;
    MeasureType m_QuantileValue;
    MeasureType m_MaximalValue;

    MeasureType m_TargetMeanValue;
    MeasureType m_TargetRMSEValue;
    MeasureType m_TargetQuantileValue;
    MeasureType m_TargetMaximalValue;

    template <typename FixedPointSetType, typename MovingPointSetType>
    void ComputeMetrics(std::vector<MeasureType> & metrics, typename FixedPointSetType::ConstPointer fixedPointSet, typename MovingPointSetType::ConstPointer movingPointSet)
    {
      typename FixedPointSetType::PointsContainer::ConstPointer fixedContainer = fixedPointSet->GetPoints();
      typename MovingPointSetType::PointsContainer::ConstPointer movingContainer = movingPointSet->GetPoints();

      typedef itk::PointsLocator<typename FixedPointSetType::PointsContainer> PointsLocatorType;
      typename PointsLocatorType::Pointer pLocator = PointsLocatorType::New();
      pLocator->SetPoints(const_cast<typename FixedPointSetType::PointsContainer*> (fixedContainer.GetPointer()));
      pLocator->Initialize();

      ListMeasureType::Pointer measures = ListMeasureType::New();

      MeasureType mean = itk::NumericTraits<MeasureType>::Zero;
      MeasureType rmse = itk::NumericTraits<MeasureType>::Zero;
      MeasureType maximal = itk::NumericTraits<MeasureType>::Zero;

      for (typename MovingPointSetType::PointsContainerConstIterator moving = movingContainer->Begin(); moving != movingContainer->End(); ++moving) {
        typename MovingPointSetType::PointType movingPoint = moving.Value();

        size_t idx = pLocator->FindClosestPoint(movingPoint);
        MeasureType distance = movingPoint.EuclideanDistanceTo(fixedPointSet->GetPoint(idx));

        measures->PushBack(distance);
        mean += distance;
        rmse += distance * distance;
        maximal = std::max(maximal, distance);
      }

      size_t numberOfPoints = movingPointSet->GetNumberOfPoints();
      mean = mean / numberOfPoints;
      rmse = std::sqrt(rmse / numberOfPoints);

      // compute quantile
      typename HistogramType::SizeType size(1);
      size.Fill(m_HistogramSize);

      HistogramFilterType::Pointer sampleToHistogram = HistogramFilterType::New();
      sampleToHistogram->SetInput(measures);
      sampleToHistogram->SetAutoMinimumMaximum(true);
      sampleToHistogram->SetHistogramSize(size);
      try {
        sampleToHistogram->Update();
      }
      catch (itk::ExceptionObject& excep) {
        itkExceptionMacro(<< excep);
      }
      MeasureType quantile = sampleToHistogram->GetOutput()->Quantile(0, m_LevelOfQuantile);

      metrics.resize(0);
      metrics.push_back(mean);
      metrics.push_back(rmse);
      metrics.push_back(quantile);
      metrics.push_back(maximal);
    }

    void ComputeTargetMetrics()
    {
      if (!m_TargetPointSet) {
        return;
      }

      typename MovingPointSetType::PointsContainer::ConstPointer movingContainer = m_MovingPointSet->GetPoints();
      typename MovingPointSetType::PointsContainerConstIterator moving = movingContainer->Begin();

      typename FixedPointSetType::PointsContainer::ConstPointer targetContainer = m_TargetPointSet->GetPoints();
      typename FixedPointSetType::PointsContainerConstIterator target = targetContainer->Begin();

      ListMeasureType::Pointer measures = ListMeasureType::New();

      m_TargetMeanValue = itk::NumericTraits<MeasureType>::Zero;
      m_TargetRMSEValue = itk::NumericTraits<MeasureType>::Zero;
      m_TargetMaximalValue = itk::NumericTraits<MeasureType>::Zero;

      for ( ; moving != movingContainer->End(); ++moving, ++target) {
        typename MovingPointSetType::PointType movingPoint = moving.Value();
        MeasureType distance = movingPoint.EuclideanDistanceTo(target.Value());

        measures->PushBack(distance);

        m_TargetMeanValue += distance;
        m_TargetRMSEValue += distance * distance;
        m_TargetMaximalValue = std::max(m_TargetMaximalValue, distance);
      }

      size_t numberOfPoints = m_TargetPointSet->GetNumberOfPoints();
      m_TargetMeanValue = m_TargetMeanValue / numberOfPoints;
      m_TargetRMSEValue = std::sqrt(m_TargetRMSEValue / numberOfPoints);

      // compute quantile
      typename HistogramType::SizeType size(1);
      size.Fill(m_HistogramSize);

      HistogramFilterType::Pointer sampleToHistogram = HistogramFilterType::New();
      sampleToHistogram->SetInput(measures);
      sampleToHistogram->SetAutoMinimumMaximum(true);
      sampleToHistogram->SetHistogramSize(size);
      sampleToHistogram->Update();
      m_TargetQuantileValue = sampleToHistogram->GetOutput()->Quantile(0, m_LevelOfQuantile);
    }
  };
}
