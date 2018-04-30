#ifndef itkGMML2RigidPointSetToPointSetMetric_h
#define itkGMML2RigidPointSetToPointSetMetric_h

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
class GMML2RigidPointSetToPointSetMetric : public GMMPointSetToPointSetMetricBase < TFixedPointSet, TMovingPointSet >
{
public:
  /** Standard class typedefs. */
  typedef GMML2RigidPointSetToPointSetMetric                                  Self;
  typedef GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >  Superclass;
  typedef SmartPointer< Self >                                                Pointer;
  typedef SmartPointer< const Self >                                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GMML2RigidPointSetToPointSetMetric, GMMPointSetToPointSetMetricBase);

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

  /** Calculates the local metric value for a single point.*/
  virtual MeasureType GetLocalNeighborhoodValue(const MovingPointIterator &) const ITK_OVERRIDE;

  /** Calculates the local value/derivative for a single point.*/
  virtual bool GetLocalNeighborhoodValueAndDerivative(const MovingPointIterator &, MeasureType &, LocalDerivativeType &) const ITK_OVERRIDE;

  /** Calculates the local derivative for a single point.*/
  virtual bool GetLocalNeighborhoodDerivative(const MovingPointIterator &, LocalDerivativeType &) const ITK_OVERRIDE;

  virtual double GetNormalizingValueFactor() const ITK_OVERRIDE
  {
    return -1.0 / (this->m_FixedPointSet->GetNumberOfPoints() * this->m_MovingPointSet->GetNumberOfPoints());
  }

  virtual double GetNormalizingDerivativeFactor() const ITK_OVERRIDE
  {
    return 2.0 / (this->m_FixedPointSet->GetNumberOfPoints() * this->m_MovingPointSet->GetNumberOfPoints()) / (this->m_Scale * this->m_Scale);
  }

  virtual double GetNormalizingDerivativeScaleFactor() const ITK_OVERRIDE
  {
    return -4.0 / (this->m_FixedPointSet->GetNumberOfPoints() * this->m_MovingPointSet->GetNumberOfPoints()) / (this->m_Scale * this->m_Scale * this->m_Scale);
  }

protected:
  GMML2RigidPointSetToPointSetMetric();
  virtual ~GMML2RigidPointSetToPointSetMetric() {}

private:
  GMML2RigidPointSetToPointSetMetric(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGMML2RigidPointSetToPointSetMetric.hxx"
#endif

#endif
