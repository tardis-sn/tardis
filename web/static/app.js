(function (){
  'use strict';
  var app = angular.module('tardisApp', ['ngResource']);

  app.config(function($interpolateProvider) {
    $interpolateProvider.startSymbol('{$');
    $interpolateProvider.endSymbol('$}');
  });

  app.config(function($logProvider){
      $logProvider.debugEnabled(true);
  });  
})();