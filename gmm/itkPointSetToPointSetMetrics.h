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

    /** Type of the additional information. */
    typedef std::pair<std::string, std::string> PairType;
    typedef std::vector<PairType> InfoType;

    void SetInfo(InfoType& info)
    {
      m_Info = info;
    }

    /** Get/Set the Fixed Point Set.  */
    itkSetConstObjectMacro(FixedPointSet, FixedPointSetType);
    itkGetConstObjectMacro(FixedPointSet, FixedPointSetType);

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
    }

    /** Compute metrics. */
    void Compute()
    {
      std::vector<MeasureType> movingToFixedMetrics;
      this->ComputeMetrics(movingToFixedMetrics, m_MovingPointSet, m_FixedPointSet);

      std::vector<MeasureType> fixedToMovingMetrics;
      this->ComputeMetrics(fixedToMovingMetrics, m_FixedPointSet, m_MovingPointSet);

      m_MeanValue = 0.5 * (movingToFixedMetrics[0] + fixedToMovingMetrics[0]);
      m_RMSEValue = 0.5 * (movingToFixedMetrics[1] + fixedToMovingMetrics[1]);
      m_QuantileValue = 0.5 * (movingToFixedMetrics[2] + fixedToMovingMetrics[2]);
      m_MaximalValue = 0.5 * (movingToFixedMetrics[3] + fixedToMovingMetrics[3]);
    }


  protected:
    PointSetToPointSetMetrics() {}
    virtual ~PointSetToPointSetMetrics() {}

  private:
    PointSetToPointSetMetrics(const Self &);
    void operator=(const Self &);

    FixedPointSetConstPointer m_FixedPointSet;
    MovingPointSetConstPointer m_MovingPointSet;

    size_t m_BucketSize = 16;
    size_t m_HistogramSize = 1000;
    double m_LevelOfQuantile = 0.95;

    MeasureType m_MeanValue;
    MeasureType m_RMSEValue;
    MeasureType m_QuantileValue;
    MeasureType m_MaximalValue;
    InfoType m_Info;

    void ComputeMetrics(std::vector<MeasureType> & metrics, typename MovingPointSetType::ConstPointer movingPointSet, typename FixedPointSetType::ConstPointer fixedPointSet)
    {
      typename FixedPointSetType::PointsContainer::ConstPointer fixedContainer = fixedPointSet->GetPoints();
      typename MovingPointSetType::PointsContainer::ConstPointer movingContainer = movingPointSet->GetPoints();

      typedef itk::PointsLocator<typename FixedPointSetType::PointsContainer> PointsLocatorType;
      typename PointsLocatorType::Pointer pLocator = PointsLocatorType::New();
      pLocator->SetPoints(const_cast<typename FixedPointSetType::PointsContainer*> (fixedContainer.GetPointer()));
      pLocator->Initialize();

      typedef itk::Vector<MeasureType, 1> VectorType;
      typedef itk::Statistics::ListSample<VectorType> ListSampleType;
      ListSampleType::Pointer measures = ListSampleType::New();

      MeasureType mean = itk::NumericTraits<MeasureType>::Zero;
      MeasureType rmse = itk::NumericTraits<MeasureType>::Zero;
      MeasureType maximal = itk::NumericTraits<MeasureType>::Zero;

      for (typename MovingPointSetType::PointsContainerConstIterator moving = movingContainer->Begin(); moving != movingContainer->End(); ++moving) {
        typename FixedPointSetType::PointType movingPoint = moving.Value();

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
      typedef typename itk::Statistics::Histogram<MeasureType, itk::Statistics::DenseFrequencyContainer2> HistogramType;
      typename HistogramType::SizeType size(1);
      size.Fill(m_HistogramSize);

      typedef itk::Statistics::SampleToHistogramFilter<ListSampleType, HistogramType> SampleToHistogramFilterType;
      SampleToHistogramFilterType::Pointer sampleToHistogram = SampleToHistogramFilterType::New();
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
  };
}
