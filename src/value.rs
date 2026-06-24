/// A single decoded INFO/FORMAT scalar.
#[derive(Debug, Clone, PartialEq)]
pub enum Scalar {
    Int(i64),
    Float(f64),
    Char(char),
    Str(String),
}

/// A decoded INFO/FORMAT field value.
#[derive(Debug, Clone, PartialEq)]
pub enum FieldValue {
    Flag,
    Scalar(Scalar),
    List(Vec<Option<Scalar>>),
}

impl FieldValue {
    /// Number of list entries, or `None` for Flag/lone-scalar values.
    pub fn list_len(&self) -> Option<usize> {
        match self {
            FieldValue::List(v) => Some(v.len()),
            _ => None,
        }
    }

    pub fn ints<I: IntoIterator<Item = i64>>(xs: I) -> FieldValue {
        FieldValue::List(xs.into_iter().map(|x| Some(Scalar::Int(x))).collect())
    }

    pub fn floats<I: IntoIterator<Item = f64>>(xs: I) -> FieldValue {
        FieldValue::List(xs.into_iter().map(|x| Some(Scalar::Float(x))).collect())
    }

    pub fn strings<I, S>(xs: I) -> FieldValue
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        FieldValue::List(
            xs.into_iter()
                .map(|x| Some(Scalar::Str(x.into())))
                .collect(),
        )
    }
}
