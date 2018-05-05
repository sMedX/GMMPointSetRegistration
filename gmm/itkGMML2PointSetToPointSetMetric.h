#ifndef itkGMML2PointSetToPointSetMetric_h
#define itkGMML2PointSetToPointSetMetric_h

#include "itkGMMPointSetToPointSetMetricBase.h"

namespace itk
{
/** \class GMMRigidPointSetToPointSetMetric
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
class GMML2PointSetToPointSetMetric : public GMMPointSetToPointSetMetricBase < TFixedPointSet, TMovingPointSet >
{
public:
  /** Standard class typedefs. */
  typedef GMML2PointSetToPointSetMetric                                       Self;
  typedef GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >  Superclass;
  typedef SmartPointer< Self >                                                Pointer;
  typedef SmartPointer< const Self >                                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GMML2PointSetToPointSetMetric, GMMPointSetToPointSetMetricBase);

  /** Types transferred from the base class */
  typedef typename Superclass::MeasureType                    MeasureType;
  typedef typename Superclass::FixedPointType                 FixedPointType;
  typedef typename Superclass::FixedPointIterator             FixedPointIterator;
  typedef typename Superclass::MovingPointType                MovingPointType;
  typedef typename Superclass::MovingPointIterator            MovingPointIterator;
  typedef typename Superclass::DerivativeValueType            DerivativeValueType;
  typedef typename Superclass::LocalDerivativeType            LocalDerivativeType;
  typedef typename Superclass::FixedNeighborsIdentifierType   FixedNeighborsIdentifierType;
  typedef typename Superclass::MovingNeighborsIdentifierType  MovingNeighborsIdentifierType;
  typedef typename Superclass::FixedNeighborsIteratorType     FixedNeighborsIteratorType;
  typedef typename Superclass::MovingNeighborsIteratorType    MovingNeighborsIteratorType;

  /** Initialize the Metric by making sure that all the components are present and plugged together correctly.*/
  virtual void Initialize() throw (ExceptionObject) ITK_OVERRIDE;

protected:
  GMML2PointSetToPointSetMetric();
  virtual ~GMML2PointSetToPointSetMetric() {};

  /** Calculates the local metric value for a single point.*/
  virtual MeasureType GetLocalNeighborhoodValue(const MovingPointIterator &) const ITK_OVERRIDE;

  /** Calculates the local value/derivative for a single point.*/
  virtual bool GetLocalNeighborhoodValueAndDerivative(const MovingPointIterator &, MeasureType &, LocalDerivativeType &) const ITK_OVERRIDE;

  /** Calculates the local derivative for a single point.*/
  virtual bool GetLocalNeighborhoodDerivative(const MovingPointIterator &, LocalDerivativeType &) const ITK_OVERRIDE;

  /** Initialize to prepare for a particular iteration, generally an iteration of optimization. */
  virtual void InitializeForIteration(const ParametersType & parameters) const;

  virtual double GetNormalizingValueFactor() const ITK_OVERRIDE
  {
    return -1.0 / (this->m_FixedPointSet->GetNumberOfPoints() * this->m_MovingPointSet->GetNumberOfPoints());
  }

  virtual double GetNormalizingDerivativeFactor() const ITK_OVERRIDE
  {
    return 2.0 / (this->m_FixedPointSet->GetNumberOfPoints() * this->m_MovingPointSet->GetNumberOfPoints()) / (this->m_Scale * this->m_Scale);
  }

private:
  GMML2PointSetToPointSetMetric(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGMML2PointSetToPointSetMetric.hxx"
#endif

#endif
