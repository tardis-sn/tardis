'use strict';
var app = angular.module('tardisApp', []);

app.config(function($interpolateProvider) {
  $interpolateProvider.startSymbol('{$');
  $interpolateProvider.endSymbol('$}');
});

app.config(function($logProvider){
    $logProvider.debugEnabled(true);
});