use std::cmp;
use std::fmt::Write;

pub struct Polynomial<T> {
    coefs_: Vec<T>,
}

impl<T> Polynomial<T>
where
    T: Copy
        + std::cmp::PartialEq
        + Default
        + std::ops::Add<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::AddAssign<<T as std::ops::Mul>::Output>,
{
    pub fn new(mut coefs: Vec<T>) -> Self {
        let coefs_ = coefs.clone();
        let mut c_it = coefs_.iter().rev().peekable();
        let mut cnt = 0;

        // count number of zero elements at the end of coefficients vector
		while c_it.peek().is_some() {
            if *c_it.next().unwrap() != T::default() {
                break;
            }
            if c_it.peek().is_none() {
                break;
            }
            cnt += 1;
        }

        // trim zero coefficients from the tail (if any), except zero polynomial case
        Polynomial {
            coefs_: coefs.drain(..(coefs.len() - cnt)).collect(),
        }
    }

    pub fn interpolate_from(_points: Vec<(T, T)>) -> Self {
        todo!();
    }

    pub fn add(&self, other: &Polynomial<T>) -> Self {
        let max_len = cmp::max(self.coefs_.len(), other.coefs_.len());
        let mut res_coefs = vec![T::default(); max_len];

		for (i, coef) in res_coefs.iter_mut().enumerate() {
            let mut rhs = T::default();
            let mut lhs = T::default();

            if let Some(coef) = self.coefs_.get(i) {
                rhs = *coef;
            }
            if let Some(coef) = other.coefs_.get(i) {
                lhs = *coef;
            }

            *coef = rhs + lhs;
        }

        Polynomial::new(res_coefs)
    }

    pub fn mul(&self, other: &Polynomial<T>) -> Self {
        let mut res_coefs = vec![T::default(); self.coefs_.len() + other.coefs_.len() - 1];

        for (i, a) in self.coefs_.iter().enumerate() {
            for (j, b) in other.coefs_.iter().enumerate() {
                res_coefs[i + j] += *a * *b;
            }
        }

        Polynomial::new(res_coefs)
    }

    pub fn eval_at(&self, x: T) -> T
    where
        <T as std::ops::Mul>::Output: std::ops::Add<T>,
    {
        let mut sum = T::default();

        for coef in self.coefs_.iter().rev() {
            sum = sum * x + *coef;
        }

        sum
    }
}

impl<T: std::fmt::Display> std::fmt::Debug for Polynomial<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut poly_str = "Poly: ".to_string();
        for (pow, coef) in self.coefs_.iter().enumerate() {
            write!(poly_str, "({}x^{}) + ", coef, pow)?;
        }

        write!(f, "{}", poly_str.trim_end_matches(" + "))
    }
}

impl<T: std::cmp::PartialEq> PartialEq for Polynomial<T> {
    fn eq(&self, other: &Self) -> bool {
        if self.coefs_.len() != other.coefs_.len() {
            return false;
        }
        for i in 0..self.coefs_.len() {
            if self.coefs_[i] != other.coefs_[i] {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_work() {
        let poly = Polynomial::new(vec![0.0, 1.0, 2.0, 3.0]);
        println!("{:?}", poly);
        let poly1 = Polynomial::new(vec![0.0, 1.0, 2.0, 3.0, 0.0, 0.0]);
        println!("{:?}", poly1);
        let poly1 = Polynomial::new(vec![0.0, 1.0, 2.0, 3.0, 0.0]);
        println!("{:?}", poly1);
        let poly1 = Polynomial::new(vec![0.0, 0.0, 0.0]);
        println!("{:?}", poly1);
        let poly1 = Polynomial::new(vec![1.0]);
        println!("{:?}", poly1);
        let poly1 = Polynomial::new(vec![0.0]);
        println!("{:?}", poly1);
    }

    #[test]
    fn sum() {
        let poly1 = Polynomial::new(vec![0.0, 1.0, 2.0]);
        let poly2 = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let poly3 = poly1.add(&poly2);
        assert_eq!(poly3, Polynomial::new(vec![1.0, 3.0, 5.0]));

        let poly1 = Polynomial::new(vec![1.0]);
        let poly2 = Polynomial::new(vec![0.0, 1.0, 3.0]);
        let poly3 = poly1.add(&poly2);
        assert_eq!(poly3, Polynomial::new(vec![1.0, 1.0, 3.0]));
        let poly3 = poly2.add(&poly1);
        assert_eq!(poly3, Polynomial::new(vec![1.0, 1.0, 3.0]));
    }

    #[test]
    fn mul() {
        let poly1 = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let poly2 = Polynomial::new(vec![1.0, 2.0, 3.0]);
        let poly3 = poly1.mul(&poly2);
        assert_eq!(poly3, Polynomial::new(vec![1.0, 4.0, 10.0, 12.0, 9.0]));

        let poly1 = Polynomial::new(vec![1, 2, 3]);
        let poly2 = Polynomial::new(vec![0]);
        let poly3 = poly1.mul(&poly2);
        assert_eq!(poly3, Polynomial::new(vec![0]));

        let poly1 = Polynomial::new(vec![1, 2, 3, 4]);
        let poly2 = Polynomial::new(vec![0, 1, 2, 3]);
        let poly3 = poly1.mul(&poly2);
        assert_eq!(poly3, Polynomial::new(vec![0, 1, 4, 10, 16, 17, 12]));

        let poly1 = Polynomial::new(vec![1, 2, 3]);
        let poly2 = Polynomial::new(vec![1, 2, 3, 4, 5]);
        let poly3 = poly1.mul(&poly2);
        assert_eq!(poly3, Polynomial::new(vec![1, 4, 10, 16, 22, 22, 15]));
    }

    #[test]
    fn eval() {
        let poly = Polynomial::new(vec![1, 2, 3]);
        assert_eq!(poly.eval_at(8), 209);

        let poly = Polynomial::new(vec![1i64, 4, 10, 16, 22, 22, 15]);
        assert_eq!(poly.eval_at(29), 9389554026);
    }

    #[test]
    fn interpolate() {
        todo!();
    }
}
