var gulp = require('gulp'),
    jade = require('gulp-jade');

// Templates
gulp.task('templates', function() {
  return gulp.src('src/index.jade')
    .pipe(jade())
    .pipe(gulp.dest('.'));
});

// Default task
gulp.task('default', function() {
    gulp.start('templates');
});
 
// Watch
gulp.task('watch', function() {
 
  // Watch jade files
  gulp.watch('src/**/*.jade', ['templates']);
 
});

