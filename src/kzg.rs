use blstrs::{G1Projective, G2Projective, Scalar, pairing};
use group::{ff, Group, Curve};
use ff::Field;
use crate::{rng::CountingRng, polynomial::Polynomial};
use std::ops::{Add, Mul};

pub struct PublicParams {
    pub n: usize,
    pub gen1: G1Projective,
    pub gen2: G2Projective,
    pub g1_elements:  Vec<G1Projective>,
    pub g2_elements: Vec<G2Projective>,
}

pub struct Witness {
    pub point: G1Projective,
    pub z: Scalar,
    pub y: Scalar,
}

pub fn scalar_powers(x: Scalar, n: usize) -> Vec<Scalar> {
    let mut result = Vec::new();
    let mut current = Scalar::one();
    for _ in 0..n {
        result.push(current);
        current *= x;
    }
    result
}

pub struct KZGPolyCommitmentScheme {
    pub pp: PublicParams
}

impl KZGPolyCommitmentScheme {
    pub fn setup(n: usize) -> Self {
        let gen1 = G1Projective::generator();
        let gen2 = G2Projective::generator();
        let secret = Scalar::random(CountingRng(5));
        let secret_powers = scalar_powers(secret, n);
        let g1_elements: Vec<G1Projective> = secret_powers.iter().map(|scalar| { gen1.mul(scalar) }).collect();
        let g2_elements: Vec<G2Projective> = secret_powers.iter().map(|scalar| { gen2.mul(scalar) }).collect();
        KZGPolyCommitmentScheme {
            pp: PublicParams{ n, gen1, gen2, g1_elements, g2_elements }
        }
    }

    pub fn commit(&self, poly: &Polynomial<Scalar>) -> G1Projective {

        if poly.coeffs.len() != self.pp.n {
            panic!("different length vectors")
        }
        let g1_elements = &self.pp.g1_elements;
        let mut result: G1Projective = G1Projective::identity();
        g1_elements.iter().zip(poly.coeffs.iter()).for_each(|(g1_element, coeff)| {
            result = result.add(g1_element.mul(coeff));
        });
        result
    }

    // The verifier creates a witness to the fact that p(z) = y
    pub fn create_witness(&self, poly: &Polynomial<Scalar>, z: Scalar) -> Witness  {
        
        // Calculate y := p(z)
        let y = poly.evaluate(z);

        // Calculate the numerator of q(X)
        let mut numerator = poly.clone();

        numerator.coeffs[0] = numerator.coeffs[0] - y;

        // Calculate the denominator of q(X)
        let mut denominator = Polynomial::new_zero_poly(self.pp.n + 1);
        denominator.coeffs[0] = Scalar::zero() - z;
        denominator.coeffs[1] = Scalar::one();

        let (quotient, _) = numerator.divide(&denominator, None);

        let g1_elements = &self.pp.g1_elements;
        let mut result: G1Projective = G1Projective::identity();
        g1_elements.iter().zip(quotient.coeffs.iter()).for_each(|(g1_element, coeff)| {
            result = result.add(g1_element.mul(coeff));
        });
        Witness {
            point: result,
            z,
            y,
        }
    }

    pub fn verify(&self, witness: &Witness, commitment: &G1Projective) -> bool {
        // calculate LHS e(proof, [s-z]_2)
        let s_point_in_g2 = self.pp.g2_elements[1];
        let z_point_in_g2 = self.pp.gen2.mul(witness.z);
        let s_minus_z = (s_point_in_g2 - z_point_in_g2).to_affine();
        let lhs = pairing(&witness.point.to_affine(), &s_minus_z);

        // calculater RHS e(C-[y]_1, H)
        let y_point_in_g1 = self.pp.gen1.mul(witness.y);
        let c_minus_y = (commitment - y_point_in_g1).to_affine();
        let rhs = pairing(&c_minus_y, &self.pp.gen2.to_affine());
        rhs == lhs
    }
}

#[test]
fn test_scalar_powers_1() {
    let expected: Vec<Scalar> = vec![Scalar::one(), Scalar::one(), Scalar::one()];
    let output = scalar_powers(Scalar::one(), 3);
    output.iter().zip(expected.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_scalar_powers_2() {
    let start: Scalar = Scalar::from(2);
    let expected = vec![Scalar::one(), Scalar::from(2), Scalar::from(4)];
    let output = scalar_powers(start, 3);
    output.iter().zip(expected.iter()).for_each(|(output, expected_output)| assert!(output == expected_output));
}

#[test]
fn test_e2e() {
    let kzg_instance = KZGPolyCommitmentScheme::setup(3);
    let poly = Polynomial::new(vec![Scalar::from(0) - Scalar::from(4), Scalar::zero(), Scalar::one()]);
    let commitment = kzg_instance.commit(&poly);
    let witness = kzg_instance.create_witness(&poly, Scalar::from(2));
    let result = kzg_instance.verify(&witness, &commitment);
    assert!(result == true);
}