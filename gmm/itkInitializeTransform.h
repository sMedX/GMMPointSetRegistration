#pragma once

#include <itkObject.h>
#include <itkTranslationTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkScaleSkewVersor3DTransform.h>

namespace itk
{
  template<typename TParametersValueType>
  class InitializeTransform : public Object
  {
  public:
    /** Standard class typedefs. */
    typedef InitializeTransform       Self;
    typedef Object                    Superclass;
    typedef SmartPointer<Self>        Pointer;
    typedef SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    itkTypeMacro(InitializeTransform, Object);

    enum Transform
    {
      Translation = 0,
      Versor3D = 1,
      Similarity = 2,
      ScaleSkewVersor3D = 3
    };

    /** typedefs */
    itkStaticConstMacro(PointDimension, unsigned int, 3U);
    static_assert(PointDimension == 3U, "Invalid dimension. Dimension 3 is supported.");

    typedef typename itk::Transform<TParametersValueType, PointDimension> TransformType;
    typedef typename TransformType::InputPointType InputPointType;
    typedef typename TransformType::OutputPointType OutputPointType;
    typedef typename TransformType::OutputVectorType OutputVectorType;
    typedef itk::Array<double> ParametersType;
    typedef itk::Array<unsigned int> ModeBoundsType;
    typedef itk::Array<double> BoundsType;

    // Get transform
    itkGetObjectMacro(Transform, TransformType);

    // Set/Get type of transform
    itkSetEnumMacro(TypeOfTransform, Transform);
    itkGetEnumMacro(TypeOfTransform, Transform);
    void SetTypeOfTransform(const size_t & type) { this->SetTypeOfTransform(static_cast<Transform>(type)); }

    itkSetMacro(RotationScale, double);
    itkGetMacro(RotationScale, double);

    itkSetMacro(TranslationScale, double);
    itkGetMacro(TranslationScale, double);

    itkSetMacro(ScalingScale, double);
    itkGetMacro(ScalingScale, double);

    itkSetMacro(SkewScale, double);
    itkGetMacro(SkewScale, double);

    // Get scales and bounds
    itkGetMacro(Scales, ParametersType);
    itkGetMacro(ModeBounds, ModeBoundsType);
    itkGetMacro(LowerBounds, ParametersType);
    itkGetMacro(UpperBounds, ParametersType);

    // Set/Get moving and fixed landmarks
    itkSetMacro(MovingLandmark, InputPointType);
    itkSetMacro(FixedLandmark, OutputPointType);

    void Initialize()
    {
      m_Center = m_MovingLandmark;
      m_Translation = m_FixedLandmark - m_MovingLandmark;

      switch (m_TypeOfTransform) {
      case Transform::Translation: {
        // Translation transform
        typedef typename itk::TranslationTransform<TParametersValueType, PointDimension> TranslationTransformType;
        typename TranslationTransformType::Pointer transform = TranslationTransformType::New();
        transform->Translate(m_Translation);

        m_Transform = transform;
        this->Allocate();

        // define scales
        m_NumberOfTranslationComponents = 3;

        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
          m_ModeBounds[count] = 0;
        }

        break;
      }
      case Transform::Versor3D: {
        // VersorRigid3DTransform
        typedef itk::VersorRigid3DTransform<TParametersValueType> VersorRigid3DTransformType;
        typename VersorRigid3DTransformType::Pointer transform = VersorRigid3DTransformType::New();
        transform->SetIdentity();
        transform->SetCenter(m_Center);
        transform->SetTranslation(m_Translation);

        m_Transform = transform;
        this->Allocate();

        // define scales
        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;

        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfRotationComponents; ++i, ++count) {
          m_Scales[count] = m_RotationScale;
          m_ModeBounds[count] = 2;
          m_LowerBounds[count] = -1;
          m_UpperBounds[count] = 1;
        }

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
          m_ModeBounds[count] = 0;
        }

        break;
      }


      case Transform::Similarity:{
        // Similarity3DTransform
        typedef itk::Similarity3DTransform<TParametersValueType> Similarity3DTransformType;
        typename Similarity3DTransformType::Pointer transform = Similarity3DTransformType::New();
        transform->SetIdentity();
        transform->SetCenter(m_Center);
        transform->SetTranslation(m_Translation);

        m_Transform = transform;
        this->Allocate();

        // define scales
        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;
        m_NumberOfScalingComponents = 1;

        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfRotationComponents; ++i, ++count) {
          m_Scales[count] = m_RotationScale;
          m_ModeBounds[count] = 2;
          m_LowerBounds[count] = -1;
          m_UpperBounds[count] = 1;
        }

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
          m_ModeBounds[count] = 0;
        }

        for (size_t i = 0; i < m_NumberOfScalingComponents; ++i, ++count) {
          m_Scales[count] = m_ScalingScale;
          m_ModeBounds[count] = 0;
        }

        break;
      }

      case Transform::ScaleSkewVersor3D:{
        typedef itk::ScaleSkewVersor3DTransform<TParametersValueType> ScaleSkewVersor3DTransformType;
        typename ScaleSkewVersor3DTransformType::Pointer transform = ScaleSkewVersor3DTransformType::New();
        transform->SetIdentity();
        transform->SetCenter(m_Center);
        transform->SetTranslation(m_Translation);

        m_Transform = transform;
        this->Allocate();

        // define scales
        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;
        m_NumberOfScalingComponents = 3;
        m_NumberOfSkewComponents = 6;

        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfRotationComponents; ++i, ++count) {
          m_Scales[count] = m_RotationScale;
          m_ModeBounds[count] = 2;
          m_LowerBounds[count] = -1;
          m_UpperBounds[count] = 1;
        }

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
          m_ModeBounds[count] = 0;
        }

        for (size_t i = 0; i < m_NumberOfScalingComponents; ++i, ++count) {
          m_Scales[count] = m_ScalingScale;
          m_ModeBounds[count] = 0;
        }

        for (size_t i = 0; i < m_NumberOfSkewComponents; ++i, ++count) {
          m_Scales[count] = m_SkewScale;
          m_ModeBounds[count] = 0;
        }

        break;
      }
      default:
        itkExceptionMacro(<< "Invalid type of the input transform");
      }
    }

  protected:
    InitializeTransform() {}
    virtual ~InitializeTransform() {};

    virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE
    {
      Superclass::PrintSelf(os, indent);
      os << std::endl;

      os << "Landmarks" << std::endl;
      os << indent << "fixed  " << m_FixedLandmark << std::endl;
      os << indent << "moving " << m_MovingLandmark << std::endl;
      os << std::endl;

      os << "Transform " << m_Transform->GetTransformTypeAsString() << std::endl;
      os << indent << "Center           " << m_Center << std::endl;
      os << indent << "Translation      " << m_Translation << std::endl;
      os << indent << "Fixed parameters " << m_Transform->GetFixedParameters() << ", " << m_Transform->GetNumberOfFixedParameters() << std::endl;
      os << indent << "Parameters       " << m_Transform->GetParameters() << ", " << m_Transform->GetNumberOfParameters() << std::endl;
      os << std::endl;
      os << indent << "Scales           " << m_Scales << std::endl;
      os << indent << "Mode bounds      " << m_ModeBounds << std::endl;
      os << indent << "Lower bounds     " << m_LowerBounds << std::endl;
      os << indent << "Upper bounds     " << m_UpperBounds << std::endl;
      os << std::endl;
    }

    Transform m_TypeOfTransform = Transform::Similarity;
    typename TransformType::Pointer m_Transform = nullptr;

    InputPointType m_Center;
    OutputVectorType m_Translation;

    InputPointType m_MovingLandmark;
    OutputPointType m_FixedLandmark;

    /** Set the boundary condition for each variable, where
    * select[i] = 0 if x[i] is unbounded,
    *           = 1 if x[i] has only a lower bound,
    *           = 2 if x[i] has both lower and upper bounds, and
    *           = 3 if x[1] has only an upper bound */
    ModeBoundsType m_ModeBounds;
    BoundsType m_LowerBounds;
    BoundsType m_UpperBounds;
    ParametersType m_Scales;

    double m_TranslationScale = 1;
    double m_RotationScale = 0.1;
    double m_ScalingScale = 0.1;
    double m_SkewScale = 0.1;

    size_t m_NumberOfComponents = 0;
    size_t m_NumberOfTranslationComponents = 0;
    size_t m_NumberOfRotationComponents = 0;
    size_t m_NumberOfScalingComponents = 0;
    size_t m_NumberOfSkewComponents = 0;
    size_t m_NumberOfParameters = 0;

    void Allocate()
    {
      m_NumberOfParameters = m_Transform->GetNumberOfParameters();

      m_Scales.set_size(m_NumberOfParameters);
      m_Scales.Fill(0);

      m_ModeBounds.set_size(m_NumberOfParameters);
      m_ModeBounds.Fill(0);

      m_LowerBounds.set_size(m_NumberOfParameters);
      m_LowerBounds.Fill(0);

      m_UpperBounds.set_size(m_NumberOfParameters);
      m_UpperBounds.Fill(0);
    }

  private:
    InitializeTransform(const Self &) ITK_DELETE_FUNCTION;
    void operator=(const Self &) ITK_DELETE_FUNCTION;
  };
}
