//! QA Passport generation — structured per-sample QC report + HTML dashboard

pub mod html;
mod passport;

pub use html::generate_html_report;
pub use passport::Passport;
