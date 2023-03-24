
// Storing Characters
// X: Timeline
// Y: X coordinates of character
// Z: Y coordinates of character
//
// ......
// ......
// ......
// ......
// ......
//
// Input:
// ..##@.
// .#....
// ..##..
// ....#.
// .###..
//
// ...1x0
// ...000
// ......
// ......
// ......
//
// ...x1.
// ......
// ......
// ......
// ......
//

use std::rc::Rc;
use image::{ImageBuffer, Luma};
use std::path::Path;
// type &Cell = Rc<RefCell<Cell>>;
type CellRef = Rc<Cell>;
type OutContainer = (CellRef, bool);

pub struct Canvas {
    grid: Grid,
    outline: Vec<CellRef>,
    skeleton: Vec<CellRef>,
}


impl Canvas {
    pub fn new() -> Canvas {
        Canvas { grid: Grid::new(), outline: Vec::new(), skeleton: Vec::new() }
    }

    // pub fn draw_from_str(&mut self, str_grid: &str) {
    //     self.grid.set_grid(str_grid);
    //     self.draw_skeleton();
    // }

    pub fn build_from_img(&mut self, str_src_path: &str) {
        println!("Reading image...");

        let src_path = Path::new(str_src_path);

        let img = image::open(src_path).expect("Error opening input image").to_luma8();

        let (width, height) = img.dimensions();
        self.grid.set_dim((width, height));
        println!("Images size (width, height): {}, {}", width, height);

        // let pixels: Vec<u8> = img
        //     .pixels()
        //     .map(|p| if p.0[0] > 128 { 1 } else { 0 })
        //     .collect();
        
        for (x, y, pixel) in img.enumerate_pixels() {
            if (x != 0 && y != 0)
                && (x != width-1 && y != height-1)
                && pixel.0[0] < 128 {
                self.grid.set_cell(x, y, true);
            }
        }

        println!("Building outline...");
        self.build_outline();
        println!("Building skeleton... (may take a while)");
        self.build_skeleton();
    

        // let mut new_img = ImageBuffer::<Luma<u8>, Vec<u8>>::new(width, height);
        //
        // for (x, y, pixel) in new_img.enumerate_pixels_mut() {
        //     let index = (y * width + x) as usize;
        //     *pixel = Luma([pixels[index] * 255]);
        // }
        //
        // new_img
        //     .save(dest_path)
        //     .expect("Error saving output image");
    }

    pub fn save(&self) {
        let original_path = Path::new("original.png");
        let outline_path = Path::new("outline.png");
        let skeleton_path = Path::new("skeleton.png");

        let (width, height) = self.grid.dim;
        // Saving original image
        let mut new_img = ImageBuffer::<Luma<u8>, Vec<u8>>::new(width, height);

        for (x, y, pixel) in new_img.enumerate_pixels_mut() {
            *pixel = match self.grid.at(x, y) {
                Some(_) => Luma([0]),
                _ => Luma([255]),
            }
        }

        new_img
            .save(original_path)
            .expect("Error saving original image");

        // Saving outline image
        let mut new_img = ImageBuffer::<Luma<u8>, Vec<u8>>::new(width, height);

        for (x, y, pixel) in new_img.enumerate_pixels_mut() {
            let in_outline = self.outline
                .iter()
                .find(|cell| cell.x == x as i32 && cell.y == y as i32);
            *pixel = match in_outline {
                Some(_) => Luma([0]),
                _ => Luma([255]),
            }
        }

        new_img
            .save(outline_path)
            .expect("Error saving outline image");

        // Saving skeleton image
        let mut new_img = ImageBuffer::<Luma<u8>, Vec<u8>>::new(width, height);

        for (x, y, pixel) in new_img.enumerate_pixels_mut() {
            let in_skeleton = self.skeleton
                .iter()
                .find(|cell| cell.x == x as i32 && cell.y == y as i32);
            *pixel = match in_skeleton {
                Some(_) => Luma([0]),
                _ => Luma([255]),
            }
        }

        new_img
            .save(skeleton_path)
            .expect("Error saving skeleton image");
    }

    fn build_outline(&mut self) {
        self.outline = self.grid.select(|cell| self.grid.get_adjacent(&cell).len() < 4);
    }

    fn build_skeleton(&mut self) {
        let mut skeleton: Vec<CellRef> = Vec::new();

        let mut outline: Vec<OutContainer> = self.outline
            .iter()
            .map(|cell| (Rc::clone(cell), true))
            .collect();
        let taken_cells = self.grid.select(|_| true);
        // Why tf do I need to clone it twice qwq
        let mut interior: Vec<CellRef> = taken_cells
            .iter()
            .map(|cell| Rc::clone(cell))
            .filter(|cell| !outline.contains(&(Rc::clone(cell), true)))
            .collect();

        let mut in_cell_dists: Vec<(CellRef, f64)> = Vec::new();

        // Iterating over all interiors to compare with the outline
        for in_cell in interior.iter() {
            // Get all dists from in_cell with every outline
            let dists: Vec<f64> = outline
                .iter()
                .map(|(out_cell, _)| out_cell.dist(Rc::clone(&in_cell)))
                .collect();


            // Finding the min dist
            let mut min_dist = dists[0];
            for &dist in dists.iter() {
                if dist < min_dist {
                    min_dist = dist;
                }
            }
            println!("{}", min_dist);
            in_cell_dists.push((Rc::clone(&in_cell), min_dist));
        }

        while interior.len() > 0 {
            // Finding the best interior cell that is farthest from all outline edges
            // (best middle point)
            //
            // AAAA- please actually find the solution to this temp fix
            let mut farthest_in_cell = (Rc::new(Cell { x: 0, y: 0}), 0f64);
            for (in_cell, dist) in in_cell_dists.iter() {
                if dist > &farthest_in_cell.1 {
                    farthest_in_cell = (Rc::clone(in_cell), *dist);
                }
            }

            let (farthest_in_cell, min_dist) = farthest_in_cell;

            // Get all dists from in_cell with every outline (again)
            // Optimization Idea: Don't do this...
            let dists: Vec<f64> = outline
                .iter()
                .map(|(out_cell, _)| out_cell.dist(Rc::clone(&farthest_in_cell)))
                .collect();

            
            // Getting the outline cells that are the min dists or 1 away from min dist
            let mut closest_out_cells: Vec<CellRef> = Vec::new();
            for (&dist, (out_cell, available)) in dists.iter().zip(outline.iter()) {
                if dist / 2f64.sqrt() < min_dist && *available {
                    println!("{}, {}", dist, min_dist);
                    closest_out_cells.push(Rc::clone(&out_cell));
                }
            }

            // Getting a pair of outline cells that are somewhat the same slope
            //
            // You need to do this because we only want pairs where you can draw a line that
            // intersects all three points
            let pair = self.filter_same_slope(Rc::clone(&farthest_in_cell), closest_out_cells);
            // Consuming the availablity on the outline cells so we cant reuse them
            // Inefficent but solving it the right way takes time that i don't have
            if let Some((out_cell_a, out_cell_b)) = &pair {
                for (cell, available) in outline.iter_mut() {
                    if cell == out_cell_a || cell == out_cell_b {
                        *available = false;
                    }
                }

                // Find and remove element to not be reused
                if let Some(index) = interior.iter().position(|in_cell| *in_cell == farthest_in_cell) {
                    skeleton.push(farthest_in_cell);
                    interior.remove(index);
                    in_cell_dists.remove(index);
                    continue;
                }
            }
            // Find and remove element to not be reused
            if let Some(index) = interior.iter().position(|in_cell| *in_cell == farthest_in_cell) {
                // Remove interior and in_cell_dists because in_ceLl_dists is being accessed
                // inside this loop
                interior.remove(index);
                in_cell_dists.remove(index);
            }
        }

        self.skeleton = skeleton;
    }

    pub fn filter_same_slope(&self, cell: CellRef, compare_cells: Vec<CellRef>) -> Option<(CellRef, CellRef)> {
        for compare_cell_a in &compare_cells {
            for compare_cell_b in &compare_cells {
                if compare_cell_a == compare_cell_b { continue }
                // println!("Cell A: {}, {}", compare_cell_a.x, compare_cell_a.y);
                // println!("Cell B: {}, {}", compare_cell_b.x, compare_cell_b.y);

                let x_len_b = (compare_cell_b.x - cell.x).abs();
                let y_len_b = (compare_cell_b.y - cell.y).abs();

                for x_offset in 0..1 {
                    for y_offset in 0..1 {
                        let x_len_a = (compare_cell_a.x - cell.x + x_offset).abs();
                        let y_len_a = (compare_cell_a.y - cell.y + y_offset).abs();

                        // println!("Offset: {}, {}", x_offset, y_offset);
                        // println!("X Lengths: {}, {}", x_len_a, x_len_b);
                        // println!("Y Lengths: {}, {}", y_len_a, y_len_b);
                        if x_len_a == 0 && x_len_b == 0 {
                            return Some((Rc::clone(compare_cell_a), Rc::clone(compare_cell_b)));
                        } else if x_len_a == 0 || x_len_b == 0 { continue }

                        // 1st grade math :p
                        let slope_a = y_len_a as f64 / x_len_a as f64;
                        let slope_b = y_len_b as f64 / x_len_b as f64;

                        if slope_a == slope_b {
                            // println!("{} == {}", slope_a, slope_b);
                            return Some((Rc::clone(compare_cell_a), Rc::clone(compare_cell_b)));
                        }
                    }
                }
            }
        }
        None
    }

    pub fn print(&self) {
        self.grid.print();
    }

    pub fn print_outline(&self) {
        self.grid.print_vertices(&self.outline);
    }

    pub fn print_skeleton(&self) {
        self.grid.print_vertices(&self.skeleton);
    }
}

pub struct Grid {
    dim: (u32, u32),
    cells: Vec<Vec<Option<CellRef>>>,
}
 
impl Grid {
    pub fn new() -> Grid {
        Grid { dim: (0, 0), cells: Vec::new() }
    }

    pub fn set_dim(&mut self, dim: (u32, u32)) {
        self.dim = dim;
        self.cells = vec![vec![None; dim.1 as usize]; dim.0 as usize];
    }

    pub fn set_grid(&mut self, str_grid: &str) {
        let str_grid_fmt: String = str_grid.chars().filter(|c| !c.is_whitespace()).collect();
        let (height, width) = self.dim;

        for row in 0..height {
            let mut row_vec: Vec<Option<CellRef>> = Vec::new();
            for col in 0..width {
                let cell = match str_grid_fmt.chars().nth((row * width + col) as usize).unwrap() {
                    // '#' => Some(Rc::new(RefCell::new(Cell { x: row as i32, y: col as i32 }))),
                    '#' => Some(Rc::new(Cell { x: row as i32, y: col as i32 })),
                    '.' => None,
                    c => panic!("Invalid symbol {}", c),
                };
                row_vec.push(cell);
            }
            self.cells.push(row_vec);
        }
    }

    pub fn set_cell(&mut self, row: u32, col: u32, some: bool) {
        let cell_opt = if some {
            Some(Rc::new(Cell { x: row as i32, y: col as i32 }))
        } else {
            None
        };

        self.cells[row as usize][col as usize] = cell_opt;
    }

    pub fn select<F>(&self, f: F) -> Vec<CellRef> where
        F: Fn(&CellRef) -> bool {

        let mut selected_cells: Vec<CellRef> = Vec::new();
        let (height, width) = self.dim;

        for row in 0..height {
            for col in 0..width {
                if let Some(taken_cell) = self.at(row as u32, col as u32) {
                    if f(&taken_cell) {
                        selected_cells.push(taken_cell);
                    }
                }
            }
        }

        selected_cells
    }

    pub fn print(&self) {
        let (height, width) = self.dim;

        for row in 0..height {
            for col in 0..width {
                print!("{}",
                    match self.at(row as u32, col as u32) {
                        Some(_) => '#',
                        None => '.'
                    }
                )
            }
            print!("\n");
        }
    }

    pub fn print_vertices(&self, vertices: &Vec<CellRef>) {
        let (height, width) = self.dim;

        let find = |x: i32, y: i32| {
            for v in vertices {
                if v.get_pos() == (x, y) {
                    return true;
                }
            }
            false
        };

        for row in 0..height {
            for col in 0..width {
                if find(row as i32, col as i32) {
                    print!("#");
                } else {
                    print!(".");
                }
            }
            print!("\n");
        }
    }

    pub fn at(&self, row: u32, col: u32) -> Option<CellRef> {
        // if let Some(row_vec) = self.cells.get(row) {
        //     if let Some(cell_opt) = row_vec.get(col) {
        //         // return Some(cell_opt.as_ref().map(Rc::clone));
        //         return Rc::new(cell_opt);
        //     }
        // }
        // None
        let cell = &self.cells[row as usize][col as usize];
        match cell {
            Some(taken) => Some(Rc::clone(taken)),
            _ => None,
        }

    }

    pub fn get_adjacent(&self, cell: &Cell) -> Vec<CellRef> {
        let mut adjacent: Vec<CellRef> = Vec::new();
        let (x, y) = cell.get_pos();

        for row_offset in -1..=1 {
            for col_offset in -1..=1 {
                // Adjecent not including the element itself, and diaginal
                if (row_offset != 0 && col_offset != 0)
                    || (row_offset == 0 && col_offset == 0) {
                    continue;
                }

                let row = x + row_offset;
                let col = y + col_offset;

                // Out of bounds
                if row < 0 || col < 0 {
                    continue;
                }

                let (row, col) = (row as usize, col as usize);

                // Out of bounds
                if row as u32 >= self.dim.0 || col as u32 >= self.dim.1 {
                    continue;
                }

                if let Some(taken_cell) = self.at(row as u32, col as u32) {
                    adjacent.push(taken_cell);
                }
            }
        }

        adjacent
    }
}

#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub struct Cell {
    x: i32,
    y: i32,
}

impl Cell {
    pub fn get_pos(&self) -> (i32, i32) {
        (self.x, self.y)
    }

    pub fn dist(&self, other: CellRef) -> f64 {
        let x_len = (self.x - other.x).abs();
        let y_len = (self.y - other.y).abs();
        ((x_len.pow(2) + y_len.pow(2)) as f64).sqrt()
    }
}


#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn same_slope() {
        let canvas = Canvas::new();

        let true_cell = Rc::new(Cell { x: 2, y: 2 });
        let true_compare_cells = vec![Rc::new(Cell { x: 0, y: 0 }), Rc::new(Cell { x: 4, y: 4 })];
        // println!("{}", canvas.filter_same_slope(true_cell, true_compare_cells).is_some());
        assert!(canvas.filter_same_slope(true_cell, true_compare_cells).is_some());

        let false_cell = Rc::new(Cell { x: 0, y: 0 });
        let false_compare_cells = vec![Rc::new(Cell { x: 4, y: 0 }), Rc::new(Cell { x: 0, y: 4 }), Rc::new(Cell { x: 4, y: 4 })];
        assert!(canvas.filter_same_slope(false_cell, false_compare_cells).is_none());
    }
}
