use std::iter::{IntoIterator, Iterator};

pub struct SeqWindowSlider<'a> {
    seq: &'a str,
    window_size: usize,
    offset: usize,
}

impl<'a> SeqWindowSlider<'a> {
    pub fn new(seq: &'a str, window_size: usize) -> Self {
        Self {
            seq,
            window_size,
            offset: 0,
        }
    }
}

pub struct SeqWindowSliderIterator<'a> {
    slider: SeqWindowSlider<'a>,
}

impl<'a> IntoIterator for SeqWindowSlider<'a> {
    type Item = &'a str;
    type IntoIter = SeqWindowSliderIterator<'a>;
    fn into_iter(self) -> Self::IntoIter {
        SeqWindowSliderIterator { slider: self }
    }
}

impl<'a> Iterator for SeqWindowSliderIterator<'a> {
    type Item = &'a str;

    fn next(&mut self) -> Option<Self::Item> {
        let slider = &mut self.slider;
        if slider.seq.len() - slider.window_size < slider.offset {
            None // iterator exhausted
        } else {
            let result = &slider.seq[slider.offset..slider.offset + slider.window_size];
            slider.offset += 1;
            Some(result)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seq_window_slider() {
        let seq = "01234567890";
        let slider = SeqWindowSlider::new(seq, 3);
        let mut seen_all = false;
        for (i, win) in slider.into_iter().enumerate() {
            match i {
                0 => assert_eq!(win, "012"),
                1 => assert_eq!(win, "123"),
                2 => assert_eq!(win, "234"),
                3 => assert_eq!(win, "345"),
                4 => assert_eq!(win, "456"),
                5 => assert_eq!(win, "567"),
                6 => assert_eq!(win, "678"),
                7 => assert_eq!(win, "789"),
                8 => {
                    assert_eq!(win, "890");
                    seen_all = true
                }
                _ => assert!(false),
            }
        }
        assert!(seen_all);
    }
}
