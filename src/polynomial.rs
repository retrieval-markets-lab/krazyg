use blstrs::Scalar;
use group::ff::Field;

use crate::kzg::scalar_powers;

#[derive(Debug)]
#[derive(Clone)]
pub struct Polynomial<T> {
    pub coeffs: Vec<T>
}


impl Polynomial<Scalar> {
    pub fn new(coeffs: Vec<Scalar>) -> Self {
        Polynomial { coeffs }
    }

    pub fn new_zero_poly(length: usize) -> Self {
        Polynomial { 
            coeffs: vec![Scalar::zero(); length],
        }
    }

    pub fn degree(&self) -> Option<usize> {
        let coeffs = &self.coeffs;
        let condition = |x: &Scalar| x.ne(&Scalar::zero());
        for (index, value) in coeffs.into_iter().enumerate().rev() {
            if condition(&value) {
                return Some(index)
            }
        }
        None
    }

    pub fn multiply(&self, other_poly: &Self) -> Self {
        // TODO case where self.coeffs.len is not longer than the sum of the degrees of the terms

        let mut result = vec![Scalar::zero(); self.coeffs.len()];
    
        for i in 0..self.coeffs.len() - 1 {
            for j in 0..other_poly.coeffs.len() - 1 {
                if self.coeffs[i] != Scalar::zero() && other_poly.coeffs[j] != Scalar::zero() {
                    result[i + j] += self.coeffs[i] * other_poly.coeffs[j];
                }
            }
        }
        Polynomial { 
            coeffs: result,
        }
    }

    pub fn subtract(&self, other_poly: &Self) -> Self {
        let mut result_coeffs: Vec<Scalar> = Vec::new();

        for i in 0..std::cmp::max(self.coeffs.len(), other_poly.coeffs.len()) {
            let x = if i < self.coeffs.len() { self.coeffs[i] } else { Scalar::zero() };
            let y = if i < other_poly.coeffs.len() { other_poly.coeffs[i] } else { Scalar::zero() };
            result_coeffs.push(x - y);
        }
        Polynomial { 
            coeffs: result_coeffs,
        }
    }

    // Divide the polynomial by another polynomial using long division
    pub fn divide(&self, divisor: &Polynomial<Scalar>, quotient: Option<Polynomial<Scalar>>) -> (Polynomial<Scalar>, Polynomial<Scalar>) {
        let zero_quotient = Polynomial::new_zero_poly(self.coeffs.len());
        let zero_remainder = Polynomial::new_zero_poly(self.coeffs.len());
        let dividend = self.clone();

        // TODO edge case where divisor is zero

        // Check if dividend is zero
        if dividend.coeffs.eq(&zero_quotient.coeffs) {
            match quotient {
                Some(poly)=> return (poly, zero_remainder),
                None => return (zero_quotient, zero_remainder),
            } 
        }
        
        // definitely succeeds following above checks;
        let degree = self.degree().unwrap();
        let divisor_degree = divisor.degree().unwrap();

        if degree < divisor_degree {
            match quotient {
                Some(poly)=> return (poly, dividend),
                None => return (zero_quotient, dividend),
            }
        }

        // 1. Divide the first term of the dividend by the highest term of the divisor
        let mut new_quotient = match quotient {
            Some(poly) => poly,
            None => zero_quotient,
        };
        let new_quotient_value = dividend.coeffs[degree] * divisor.coeffs[divisor_degree].invert().unwrap();
        new_quotient.coeffs[degree - 1] = new_quotient_value;
        let mut new_multiplier = Polynomial::new_zero_poly(self.coeffs.len());
        new_multiplier.coeffs[degree - 1] = new_quotient_value;

        // 2. Multiply the divisor by the result just obtained
        let multiple = divisor.multiply(&new_multiplier);

        // 3. Subtract the product just obtained from the appropriate terms of the original dividend
        let remainder = dividend.subtract(&multiple);

        // 4. Repeat the previous three steps
        remainder.divide(divisor, Some(new_quotient))
    }

    pub fn evaluate(&self, value: Scalar) -> Scalar {
        let mut result = Scalar::zero();
        let powers_of_z = scalar_powers(value, self.coeffs.len() + 1);
        self.coeffs.iter().zip(powers_of_z.iter()).for_each(|(coeff, power)| {
            result += coeff * power
        });
        result
    }
}



