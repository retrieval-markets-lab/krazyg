use blstrs::Scalar;
use group::ff::Field;

use crate::kzg::scalar_powers;

#[derive(Debug)]
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

    pub fn clone(&self) -> Self {
        let coeffs = self.coeffs.clone();
        Polynomial {
            coeffs
        }
    }

    pub fn multiply(&self, other_poly: &Self) -> Self {
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

        if self.coeffs.eq(&zero_quotient.coeffs) {
            match quotient {
                Some(poly)=> return (poly, zero_remainder),
                None => return (zero_quotient, zero_remainder),
            } 
        }
        
        let dividend = Polynomial::new(self.clone().coeffs);

        let degree = self.degree().unwrap();

        if degree < divisor.degree().unwrap() {
            match quotient {
                Some(poly)=> return (poly, dividend),
                None => return (zero_quotient, dividend),
            }
        }

        // 1. Divide the first term of the dividend by the highest term of the divisor (meaning the one with the highest power of x, which in this case is x).
        let mut new_quotient = match quotient {
            Some(poly) => poly,
            None => zero_quotient,
        };
        let new_quotient_value = dividend.coeffs[degree];
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

#[test]
fn test_divide_with_remainder() {
    // Example from https://en.wikipedia.org/wiki/Polynomial_long_division
    let mut divisor = Polynomial::new_zero_poly(4);
    divisor.coeffs[0] = Scalar::from(0) - Scalar::from(3);
    divisor.coeffs[1] = Scalar::from(1);
    let poly = Polynomial::new(vec![Scalar::from(0) - Scalar::from(4), Scalar::zero(), Scalar::from(0) - Scalar::from(2), Scalar::one()]);
    let (quotient, remainder) = poly.divide(&divisor, None);
    let expected_quotient= Polynomial::new(vec![Scalar::from(3), Scalar::one(), Scalar::one(), Scalar::zero()]);
    let expected_remainder= Polynomial::new(vec![Scalar::from(5), Scalar::zero(), Scalar::zero(), Scalar::zero()]);
    quotient.coeffs.iter().zip(expected_quotient.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
    remainder.coeffs.iter().zip(expected_remainder.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_divide_with_no_remainder() {
    // Example p(x) = (x+1)(x+2)
    let mut divisor = Polynomial::new_zero_poly(3);
    divisor.coeffs[0] = Scalar::from(1);
    divisor.coeffs[1] = Scalar::from(1);
    let poly = Polynomial::new(vec![Scalar::from(2), Scalar::from(3), Scalar::one()]);
    let (quotient, remainder) = poly.divide(&divisor, None);
    let expected_quotient= Polynomial::new(vec![Scalar::from(2), Scalar::one(),Scalar::zero(), Scalar::zero()]);
    let expected_remainder= Polynomial::new(vec![Scalar::zero(), Scalar::zero(), Scalar::zero(), Scalar::zero()]);
    quotient.coeffs.iter().zip(expected_quotient.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
    remainder.coeffs.iter().zip(expected_remainder.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}
#[test]
fn test_divide_with_no_remainder_2() {
    // Example p(x) = (x+1)(x+2)
    let mut divisor = Polynomial::new_zero_poly(3);
    divisor.coeffs[0] = Scalar::from(0) - Scalar::from(2);
    divisor.coeffs[1] = Scalar::from(1);
    let poly = Polynomial::new(vec![Scalar::from(0) - Scalar::from(4), Scalar::zero(), Scalar::one()]);
    let (quotient, remainder) = poly.divide(&divisor, None);
    let expected_quotient= Polynomial::new(vec![Scalar::from(2), Scalar::one(), Scalar::zero(), Scalar::zero()]);
    let expected_remainder= Polynomial::new(vec![Scalar::zero(), Scalar::zero(), Scalar::zero(), Scalar::zero()]);
    quotient.coeffs.iter().zip(expected_quotient.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
    remainder.coeffs.iter().zip(expected_remainder.coeffs.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_find_degree() {
    let poly_1 = Polynomial::new(vec![Scalar::one(), Scalar::one(), Scalar::one()]);
    let poly_2 = Polynomial::new(vec![Scalar::one(), Scalar::one(), Scalar::zero()]);
    let poly_3 = Polynomial::new(vec![Scalar::one(), Scalar::zero(), Scalar::zero()]);
    let expected_1 = Some(2 as usize);
    let expected_2 = Some(1 as usize);
    let expected_3 = Some(0 as usize);
    let output_1 = poly_1.degree();
    let output_2 = poly_2.degree();
    let output_3 = poly_3.degree();
    assert!(output_1 == expected_1);
    assert!(output_2 == expected_2);
    assert!(output_3 == expected_3);
}
#[test]
fn test_find_degree_when_zero() {
    let zero_poly = Polynomial::new(vec![Scalar::zero(), Scalar::zero(), Scalar::zero()]);
    assert!(zero_poly.degree() == None);
}