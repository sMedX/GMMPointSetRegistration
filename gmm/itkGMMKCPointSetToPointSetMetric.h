#ifndef itkGMMKCPointSetToPointSetMetric_h
#define itkGMMKCPointSetToPointSetMetric_h

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
class GMMKCPointSetToPointSetMetric : public GMMPointSetToPointSetMetricBase < TFixedPointSet, TMovingPointSet >
{
public:
  /** Standard class typedefs. */
  typedef GMMKCPointSetToPointSetMetric                                      Self;
  typedef GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet > Superclass;
  typedef SmartPointer< Self >                                               Pointer;
  typedef SmartPointer< const Self >                                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GMMKCPointSetToPointSetMetric, GMMPointSetToPointSetMetricBase);

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

  /** Calculates the local metric value for a single point.*/
  virtual MeasureType GetLocalNeighborhoodValue(const MovingPointType & point) const ITK_OVERRIDE;

  /** Calculates the local value/derivative for a single point.*/
  virtual void GetLocalNeighborhoodValueAndDerivative(const MovingPointType &, MeasureType &, LocalDerivativeType &) const ITK_OVERRIDE;

  /** Initialize the Metric by making sure that all the components are present and plugged together correctly.*/
  virtual void Initialize() throw (ExceptionObject) ITK_OVERRIDE;

protected:
  GMMKCPointSetToPointSetMetric();
  virtual ~GMMKCPointSetToPointSetMetric() {}

private:
  GMMKCPointSetToPointSetMetric(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGMMKCPointSetToPointSetMetric.hxx"
#endif

#endif
